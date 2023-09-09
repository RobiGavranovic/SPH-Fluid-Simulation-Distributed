import mpi.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Main {
    static List<Particle> particles;
    static Grid[][] grid;
    //neighboring particles of current processes
    static List<Particle> leftBorderParticles = new ArrayList<>();
    static List<Particle> rightBorderParticles = new ArrayList<>();
    //left and right border particles combined
    static List<List<Particle>> borderingParticles;

    //neighboring particles of adjacent processes - ghost particles
    static List<Particle> receivedLeftBorderParticles = new ArrayList<>();
    static List<Particle> receivedRightBorderParticles = new ArrayList<>();

    static ArrayList<Neighbor> neighbors = new ArrayList<>();

    //particles that become out of bounds due to movement
    static List<Particle> oobLeftBorderParticles = new ArrayList<>();
    static List<Particle> oobRightBorderParticles = new ArrayList<>();

    //particles that became in bounds due to movement
    static List<Particle> ibLeftBorderParticles = new ArrayList<>();
    static List<Particle> ibRightBorderParticles = new ArrayList<>();

    public static void main(String[] args) throws Exception {
        MPI.Init(args);
        //id of root process
        int root = 0;

        //each process has it's own rank
        int rank = MPI.COMM_WORLD.Rank();

        //number of processes (max(rank+1))
        int size = MPI.COMM_WORLD.Size();

        if (size > 20) {
            System.out.println("Too many processes, MAX (MPI.COMM_WORLD.Size) = 20. This is due to window size and gridcell size");
            System.exit(0);
        }

        //set initial particle count
        int particleCount = 3000;

        //get particles per process - to distribute them at the start evenly
        int particleCountPerProcess = particleCount / size;

        //if the division above isnt whole, add the missed out particles to root process
        if (rank == root) {
            if (particleCount % size != 0) particleCountPerProcess += particleCount % size;
        }

        //print starting statement
        if (rank == root) {
            System.out.println("Running distributed model of SPH Fluid Simulation");
            System.out.print(" with " + particleCount + " particles.\n");
        }

        //get computing area for each process
        int computingChunk = Physics.width / size;
        int fromX = rank * computingChunk;
        int toX = fromX + computingChunk;

        //now that we know which area each process is going to work on, we can initialize particles for each process in it's area
        particles = initializeParticles(particleCountPerProcess, fromX, toX);

        //initialize empty starting Grid
        grid = updateGridSize(computingChunk);

        int count = 0;
        while (count < 3000) {
            count++;

            //override lists, let garbage collector work
            oobLeftBorderParticles = new ArrayList<>();
            oobRightBorderParticles = new ArrayList<>();

            //clear neighbors
            neighbors.clear();

            //add particles to the grid
            updateGrids(fromX, toX);

            //get data of each process's particles that will have to be shared (neighboring particles)
            getMyNeighboringParticles(rank, size);

            //data of bordering collum grid's particles
            byte[] serializedRight = serialize(rightBorderParticles);
            byte[] serializedLeft = serialize(leftBorderParticles);

            //serialized "bordering collum grid's particles" length - used as information for buffer size (*particleLengthInBytes)
            byte[] serializedRightLength = serialize(serializedRight.length);
            byte[] serializedLeftLength = serialize(serializedLeft.length);

            //particle size received from neighboring processes
            int receivedParticleSizeRight = 0;
            int receivedParticleSizeLeft = 0;

            //first send information about how much data we're actually going to send

            //send right, receive from left
            if (rank == 0)
                sendParticles(serializedRightLength, rank + 1);
            else if (rank == size - 1) receivedParticleSizeLeft = receiveBufferSize(rank - 1);
            else receivedParticleSizeLeft = sendRecvBufferSize(serializedRightLength, rank + 1, rank - 1);

            //send left, receive from right
            if (rank == 0)
                receivedParticleSizeRight = receiveBufferSize(rank + 1);
            else if (rank == size - 1)
                sendParticles(serializedLeftLength, rank - 1);
            else receivedParticleSizeRight = sendRecvBufferSize(serializedLeftLength, rank - 1, rank + 1);

            //get buffer size for receiving actual data from both borders
            int bufferSizeRight = getBufferSize(receivedParticleSizeRight);
            int bufferSizeLeft = getBufferSize(receivedParticleSizeLeft);

            //each process sends rightBorderParticles to the neighboring process on the right, and receives from the process on the left.
            if (rank == 0)
                sendParticles(serializedRight, rank + 1);
            else if (rank == size - 1)
                receivedLeftBorderParticles = receiveParticles(rank - 1, bufferSizeLeft);
            else
                receivedLeftBorderParticles = sendRecvParticles(serializedRight, rank + 1, rank - 1, bufferSizeLeft);

            //each process sends leftBorderParticles to the neighboring process on the left, and receives from the process on the right.
            if (rank == 0)
                receivedRightBorderParticles = receiveParticles(rank + 1, bufferSizeRight);
            else if (rank == size - 1)
                sendParticles(serializedLeft, rank - 1);
            else
                receivedRightBorderParticles = sendRecvParticles(serializedLeft, rank - 1, rank + 1, bufferSizeRight);

            //find neighbors in grid, including the ghost particles
            findNeighbors();

            //calculate pressure for each particle
            calculatePressure();

            //calculate force for each particle
            calculateForce();

            for (Particle particle : particles) particle.move();

            //now check if particles have moved out of the process's computing area - if so then send it to another process. Receive particles from other processes

            //check out of bounds particles
            for (Particle particle : particles) {
                if (particle.x < fromX) {
                    oobLeftBorderParticles.add(particle);
                } else if (particle.x > toX) {
                    oobRightBorderParticles.add(particle);
                }
            }

            //remove out of bounds particles
            particles.removeAll(oobLeftBorderParticles);
            particles.removeAll(oobRightBorderParticles);

            //serialized particle lists to be sent off to another process
            byte[] serializedOobRight = serialize(oobRightBorderParticles);
            byte[] serializedOobLeft = serialize(oobLeftBorderParticles);

            //serialized length of data
            byte[] serializedOobRightLength = serialize(serializedOobRight.length);
            byte[] serializedOobLeftLength = serialize(serializedOobLeft.length);

            //first send information about how much data we're actually going to send
            if (rank == 0)
                sendParticles(serializedOobRightLength, rank + 1);
            else if (rank == size - 1) receivedParticleSizeLeft = receiveBufferSize(rank - 1);
            else receivedParticleSizeLeft = sendRecvBufferSize(serializedOobRightLength, rank + 1, rank - 1);

            if (rank == 0)
                receivedParticleSizeRight = receiveBufferSize(rank + 1);
            else if (rank == size - 1)
                sendParticles(serializedOobLeftLength, rank - 1);
            else receivedParticleSizeRight = sendRecvBufferSize(serializedOobLeftLength, rank - 1, rank + 1);

            //get buffer size of incoming data
            bufferSizeRight = getBufferSize(receivedParticleSizeRight);
            bufferSizeLeft = getBufferSize(receivedParticleSizeLeft);

            //send to right, receive from left
            if (rank == 0)
                sendParticles(serializedOobRight, rank + 1);
            else if (rank == size - 1) ibLeftBorderParticles = receiveParticles(rank - 1, bufferSizeLeft);
            else ibLeftBorderParticles = sendRecvParticles(serializedOobRight, rank + 1, rank - 1, bufferSizeLeft);

            //send to left, receive from right
            if (rank == 0)
                ibRightBorderParticles = receiveParticles(rank + 1, bufferSizeRight);
            else if (rank == size - 1)
                sendParticles(serializedOobLeft, rank - 1);
            else ibRightBorderParticles = sendRecvParticles(serializedOobLeft, rank - 1, rank + 1, bufferSizeRight);

            //add new in bounds particles to the list of particles
            particles.addAll(ibLeftBorderParticles);
            particles.addAll(ibRightBorderParticles);

            //if (rank == 0) System.out.println("PARTICLES: " + particles.size());

            MPI.COMM_WORLD.Barrier();

            }

        MPI.Finalize();
    }

    static List<Particle> initializeParticles(int n, int from, int to) {
        List<Particle> particles = new ArrayList<>();
        int increment = 5; //spacing between particles

        //so they arent right at the border
        from += 2;
        to -= 2;

        //put them a bit in the air
        int positionY = 10;
        int positionX = from;

        //go through all particles
        for (int i = 0; i < n; i++) {
            //if particle's X position is over process's computational area, X = from, increment Y
            if (positionX > to) {
                positionX = from;
                positionY += increment;
            }
            particles.add(new Particle(positionX, positionY));
            positionX += increment;
        }
        //return particle list, different for each process
        return particles;
    }

    // Grid Updater: Updates the size of the grid and initializes it.
    static Grid[][] updateGridSize(double width) {
        double gridSize = Physics.gridSize;
        int gridLengthX = (int) Math.floor(width / gridSize);
        Grid[][] grid = new Grid[(int) Math.floor(Physics.height / gridSize) + 1][gridLengthX];

        for (Grid[] gridArray : grid) {
            for (int i = 0; i < gridLengthX; i++) {
                gridArray[i] = new Grid();
            }
        }
        return grid;
    }

    //parameter from - because we dont have one grid anymore, but rather more grids, each starting from 0, we need to subtract the initial from value
    static void updateGrids(int fromX, int toX) {
        int gridSize = Physics.gridSize;
        int compArea = toX - fromX;

        for (Grid[] grids : grid) for (Grid grid : grids) grid.clearGrid();

        for (Particle particle : particles) {
            particle.forceX = particle.forceY = particle.density = 0;
            particle.gridX = (int) Math.floor((particle.x - fromX) / gridSize);
            particle.gridY = (int) Math.floor(particle.y / gridSize);
            if (particle.gridX < 0) particle.gridX = 0;
            if (particle.gridY < 0) particle.gridY = 0;
            if (particle.gridX > (compArea / gridSize) - 1) particle.gridX = (int) (compArea / gridSize) - 1;
            if (particle.gridY > (Physics.height / gridSize) - 1)
                particle.gridY = (int) (Physics.height / gridSize) - 1;
            grid[particle.gridY][particle.gridX].addParticle(particle);
        }
    }

    public static void getMyNeighboringParticles(int rank, int size) {
        // clear from previous iteration
        leftBorderParticles.clear();
        rightBorderParticles.clear();

        // if not rank 0, add particles from leftmost cells
        if (rank != 0) {
            for (Grid[] grids : grid) {
                leftBorderParticles.addAll(grids[0].getParticlesInGrid());
            }
        }

        // if not max rank, add particles from rightmost cells
        if (rank != size - 1) {
            for (Grid[] grids : grid) {
                rightBorderParticles.addAll(grids[grid[0].length - 1].getParticlesInGrid());
            }
        }
    }

    //input a recieved particle size and get buffer size to recieve
    public static int getBufferSize(int recievedParticleSize) {
        int bufferSizeRight;
        // 58 - null list, 385 - list and one particle, particle within list = 106, 279 - initialized list
        if (recievedParticleSize == 0) bufferSizeRight = 59;
        else bufferSizeRight = 280 + recievedParticleSize * 106;

        return bufferSizeRight;
    }

    //convert: particles to bytes
    public static byte[] serialize(Object obj) throws IOException {
        try (ByteArrayOutputStream bos = new ByteArrayOutputStream();
             ObjectOutput out = new ObjectOutputStream(bos)) {
            out.writeObject(obj);

            return bos.toByteArray();
        }
    }

    //convert: bytes to particles
    public static Object deserialize(byte[] bytes) throws IOException, ClassNotFoundException {
        try (ByteArrayInputStream bis = new ByteArrayInputStream(bytes);
             ObjectInput in = new ObjectInputStream(bis)) {

            return in.readObject();
        }
    }

    //send already serialized data to destRank
    public static void sendParticles(byte[] data, int destRank) {
        MPI.COMM_WORLD.Send(data, 0, data.length, MPI.BYTE, destRank, 0);
    }

    //receive particles list from sourceRank
    public static List<Particle> receiveParticles(int sourceRank, int bufferSize) throws IOException, ClassNotFoundException {
        byte[] buffer = new byte[bufferSize];
        MPI.COMM_WORLD.Recv(buffer, 0, buffer.length, MPI.BYTE, sourceRank, 0);

        return (List<Particle>) deserialize(buffer);
    }

    //receive information about the size of the buffer
    public static int receiveBufferSize(int sourceRank) throws IOException, ClassNotFoundException {
        byte[] buffer = new byte[100];
        MPI.COMM_WORLD.Recv(buffer, 0, buffer.length, MPI.BYTE, sourceRank, 0);

        return (int) deserialize(buffer);
    }

    //send particles list to destRank and receive particles list from SourceRank
    public static List<Particle> sendRecvParticles(byte[] data, int destRank, int sourceRank, int bufferSize) throws IOException, ClassNotFoundException {
        byte[] recvBuffer = new byte[bufferSize];

        MPI.COMM_WORLD.Sendrecv(
                data, 0, data.length, MPI.BYTE, destRank, 0,
                recvBuffer, 0, recvBuffer.length, MPI.BYTE, sourceRank, 0
        );

        return (List<Particle>) deserialize(recvBuffer);
    }

    //send buffer size to destRank and receive buffer size from SourceRank
    public static int sendRecvBufferSize(byte[] data, int destRank, int sourceRank) throws IOException, ClassNotFoundException {
        byte[] recvBuffer = new byte[100];

        MPI.COMM_WORLD.Sendrecv(
                data, 0, data.length, MPI.BYTE, destRank, 0,
                recvBuffer, 0, recvBuffer.length, MPI.BYTE, sourceRank, 0
        );

        return (int) deserialize(recvBuffer);
    }

    //findNeighbors: Clears previous neighbors list and finds current neighbors for each particle in the simulation
    static void findNeighbors() {
        neighbors.clear();
        int gridSize = Physics.gridSize;

        for (Particle particle : particles) {
            int gridX = particle.gridX;
            int gridY = particle.gridY;

            findNeighborsInGrid(particle, grid[gridY][gridX]);
            try {
                int maxX = (int) (Physics.width / gridSize) - 1;
                int maxY = (int) (Physics.height / gridSize) - 1;

                if (gridX == 0) findNeighborsInGhostParticles(particle, receivedLeftBorderParticles);
                if (gridX == maxX) findNeighborsInGhostParticles(particle, receivedRightBorderParticles);
                if (gridX < maxX) findNeighborsInGrid(particle, grid[gridY][gridX + 1]);
                if (gridY > 0) findNeighborsInGrid(particle, grid[gridY - 1][gridX]);
                if (gridX > 0) findNeighborsInGrid(particle, grid[gridY][gridX - 1]);
                if (gridY < maxY) findNeighborsInGrid(particle, grid[gridY + 1][gridX]);
                if (gridX > 0 && gridY > 0) findNeighborsInGrid(particle, grid[gridY - 1][gridX - 1]);
                if (gridX > 0 && gridY < maxY) findNeighborsInGrid(particle, grid[gridY + 1][gridX - 1]);
                if (gridX < maxX && gridY > 0) findNeighborsInGrid(particle, grid[gridY - 1][gridX + 1]);
                if (gridX < maxX && gridY < maxY) findNeighborsInGrid(particle, grid[gridY + 1][gridX + 1]);
            } catch (ArrayIndexOutOfBoundsException e) {
            }
        }
    }

    //findNeighborsInGrid: Finds neighboring particles within a specific grid cell for a given particle.
    static void findNeighborsInGrid(Particle particle, Grid gridCell) {
        double range = Physics.range;
        for (Particle particleA : gridCell.getParticlesInGrid()) {
            if (particle.equals(particleA)) continue;
            double distance = Math.pow(particle.x - particleA.x, 2) + Math.pow(particle.y - particleA.y, 2);
            if (distance < range * range) {
                Neighbor newNeighbor = new Neighbor();
                newNeighbor.setNeighbor(particle, particleA);
                neighbors.add(newNeighbor);
            }
        }
    }

    static void findNeighborsInGhostParticles(Particle particle, List<Particle> ghostParticles) {
        double range = Physics.range;
        for (Particle ghostParticle : ghostParticles) {
            double distance = Math.pow(particle.x - ghostParticle.x, 2) + Math.pow(particle.y - ghostParticle.y, 2);
            if (distance < range * range) {
                Neighbor newNeighbor = new Neighbor();
                newNeighbor.setNeighbor(particle, ghostParticle);
                neighbors.add(newNeighbor);
            }
        }
    }

    //calculatePressure: Calculates the pressure for each particle in the simulation based on its density.
    public static void calculatePressure() {
        double density = Physics.density;
        for (Particle particle : particles) {
            if (particle.density < density) particle.density = density;
            particle.pressure = particle.density - density;
        }
    }

    //calculateForce: Calculates the forces between particles based on their neighboring relationships.
    public static void calculateForce() {
        for (Neighbor neighbor : neighbors) neighbor.calculateForce();
    }
}
