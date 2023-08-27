import mpi.*;

import java.util.ArrayList;
import java.util.List;

public class Main {
    public static List<Particle> particles;

    public static void main(String args[]) throws Exception {
        MPI.Init(args);
        //id of root process
        int root = 0;

        //each process has it's own rank
        int rank = MPI.COMM_WORLD.Rank();
        //number of processes (max(rank-1))
        int size = MPI.COMM_WORLD.Size();

        //set initial particle count
        int particleCount = 100;

        //get particles per process - to distribute them at the start evenly
        int particleCountPerProcess = particleCount/size;

        //if the division above isnt whole, add the missed out particles to root process
        if (rank == root){
            if (particleCount % size != 0) particleCountPerProcess += particleCount % size;
        }

        //todo variables that might need to change destination


        //print start statement
        if (rank == root){
            System.out.println("Running distributed model of SPH Fluid Simulation");
            System.out.print(" with " + particleCount + " particles.\n");
        }

        //get computing area for each process
        int computingChunk = Physics.width / size;
        int fromX = rank * computingChunk;
        int toX = fromX + computingChunk;

        //now that we know which area each process is going to work on, we can initialize particles for each process in it's area
        particles = initializeParticles(particleCountPerProcess, fromX, toX);




        MPI.Finalize();
    }

    static List<Particle> initializeParticles(int n, int from, int to){
        List<Particle> particles = new ArrayList<>();
        int increment = 5;

        //so they arent right at the border)
        from += 2;
        to -= 2;

        //put them a bit in the air
        int positionY = 10;
        int positionX = from;

        //go trough all particles
        for (int i = 0; i < n; i++) {
            //if particle's X position is over process's comupational area, X = from, increment Y
            if (positionX > to){
                positionX = from;
                positionY += increment;
            }
            particles.add(new Particle(positionX, positionY));
            positionX += increment;
        }
        //return particle list, different for each process
        return particles;
    }
}
