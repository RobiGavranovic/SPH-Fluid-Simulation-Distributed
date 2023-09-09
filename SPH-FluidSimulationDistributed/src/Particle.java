import java.io.Serial;
import java.io.Serializable;
import java.util.Random;

public class Particle implements Serializable {
    @Serial
    private static final long serialVersionUID = 1L;

    public double velocityX, velocityY;
    public double x, y;
    public int radius;
    public int gridX, gridY;
    public double forceX, forceY;
    public double pressure, density;

    Random random = new Random();

    //todo change this variable destination - needed for move()

    public Particle(double x, double y) {
        this.x = x;
        this.y = y;
        radius = 1;

        //initial velocity so they arent falling straight down - unusual
        velocityX = random.nextInt(1 + 1) - 1;
        velocityY = random.nextInt(1 + 1) - 1;
    }

    //credit: Mitchell Sayer
    public void move(){
        this.velocityY -= Physics.gravity;
        this.velocityX += this.forceX;
        this.velocityY += this.forceY;
        this.x += this.velocityX;
        this.y += this.velocityY;

        if (this.x < 0)
            this.velocityX += (radius - this.x) * 0.5 - this.velocityX * 0.5;
        if (this.y < 0)
            this.velocityY += (radius - this.y) * 0.5 - this.velocityY * 0.5;
        if (this.x >= Physics.width)
            this.velocityX += (Physics.width - this.x) * 0.5 - this.velocityX * 0.5;
        if (this.y >= Physics.height)
            this.velocityY += (Physics.height - this.y) * 0.5 - this.velocityY * 0.5;

        //cant lets them get out of window bounds, there are no process computing areas, thus there was an issue of loss of particles
        if(this.x < 0) this.x = 1;
        if (this.x > Physics.width) x = Physics.width -1;
        if (this.y < 0) this.y = 1;
        if (this.y > Physics.height) this.y = Physics.height - 1;

    }
}
