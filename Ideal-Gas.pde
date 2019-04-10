// Make it so there are collisions until an equilibrium is reached
//(0nCollisions as the measure maybe?)

// then turn off collisions - maybe then turn on stats or have them on
// all the time to see distribution moving

int N = 50000;
int M = 1000;

float ratio = float(M)/float(N);
int setFrameRate = 60;
float dRad = 10;
float dDensity = 10;
int randVelAdd = 1;
float maxVel = 2;
float depth = 500;

float maxSpSq = 4.5;

float maxPressure = 0;
float energyMax = 0;
int maxCollisions = 0;
float t = 0;
float piSpeedSq;
float tSpeed = 0;
int nCollisions = 0;
float energy;

int count = 0;
int averagingTime = 1;


float rotAngleX, rotAngleY = 0;
PVector zero = new PVector(0, 0, 0);
float force = 0;
float pressure = 0;

float x, y, z, vx, vy, vz, k;
float posDiffSq;
//int size;
Particle pi, pj;
ArrayList<Particle> p = new ArrayList<Particle>();
int[] v = new int [M];
PVector posDiff;

void setup()
{  
  size(500, 500);//, P3D);
  frameRate(setFrameRate);

  //background(0);

  for (int i = 0; i < M; i ++)
  {
    v[i] = 0;
  }
  //v = null;

  // Make N particles with random positions and velocities
  // Add them to the ArrayList p
  for (int i = 0; i < N; i = i + 1) 
  {
    x = random(width);
    y = random(height);
    z = random(depth);

    vx = random(-randVelAdd, randVelAdd);
    vy = random(-randVelAdd, randVelAdd);
    vz = random(-randVelAdd, randVelAdd);

    //k = random(0.5, 2);

    p.add(new Particle(x, y, z, vx, vy, vz, dRad, dDensity));
    energy += x*x + y*y + z*z;
    //    size = p.size();
  }
  energy /= N;
}



void draw()
{ 
  background(0);
  //lights();


  //rotation();
  //borders();


  // set running total of force to 0 at the start of each loop

  if (count == 0)
  {
    force = 0;
    nCollisions = 0;
    tSpeed = 0;
  }
  count =+ 1;

  // Loop is how many times the system evolves per draw loop
  for (int loop = 1; loop > 0; loop = loop - 1)
  {
    // Looping through the ArrayList
    for (int i = N-1; i >= 0; i = i - 1) 
    {

      // looping through the remaining particles for each particle
      for (int j = N-1; j >= 0; j = j - 1) 
      {
        if (i != j)
        {
          // shortening the object name for the use in the loop
          pi = p.get(i);
          pj = p.get(j);

          // relative vector and constant values for the displacements
          posDiff = pi.pos.copy().sub(pj.pos.copy());
          posDiffSq = magSq(posDiff);

          // Check to see if particles are close enough to interact mechanically
          if (posDiffSq <= 4*dRad*dRad)//(pi.rad + pj.rad)*(pi.rad + pj.rad))
          {
            collide(pi, pj);
            nCollisions ++;
          } else
          {
            // Apply forces

            //gravity(pi, pj);
            //Electric(pi, pj);
          }
        }
      }
    }

    // Different loop to not have conflicts while doing calculations


    for (int i = 0; i < N; i = i + 1) 
    {
      Particle particle = p.get(i);

      // allow infinite space by wrap-around
      particle.edges();

      // apply all forces and changes
      particle.update();

      //p.get(i).show();
    }
  }

  background(0);
  if (count == averagingTime)
  {
    count = 0;

    // Calculating the pressure from the wall collisions
    // Could save this data and plot it - Thermodynamic limit
    pressure = force;//(6*width*height)*1000;
    energy = tSpeed;

    // Show all the particles in their new positons 
    for (int i = 0; i < N; i = i + 1)
    {
      //p.get(i).show();
      piSpeedSq = magSq(p.get(i).vel);
      tSpeed += piSpeedSq;

      int index = floor(M * piSpeedSq / maxSpSq);
      //println(index);
      v[index] += 1;
    }

    if (pressure > maxPressure)
    {
      maxPressure = pressure;
    }
    if (energy > energyMax)
    {
      energyMax = energy;
    }
    if (nCollisions > maxCollisions)
    {
      maxCollisions = nCollisions;
    }


    //stroke(255, 0, 0); // Red Pressure
    //point(t, height*(1-pressure/maxPressure));

    //stroke(0, 255, 0); // Green Energy
    //point(t, height*(1-energy/energyMax));

    //stroke(0, 0, 255); // Blue Collisions
    //point(t, height*(1-float(nCollisions)/float(maxCollisions)));

    //println("FrameRate: " + frameRate);
    //println("Pressure: " + pressure);
    println("Collisions: " + nCollisions);
    //println("Average Speed Squared: " + energy);
    //println("N: " + N);
    //println();

    stroke(255);
    for (int i = 0; i < M; i ++)
    {
      if (v[i] != 0)
      {
      println("count " + i + ": " + v[i]);
      }
      point(i*width/M, height*(1 - 0.1*ratio*(v[i])));
    }

    if (t > width)
    {
      t = 0;
    }
    t += 1;

    for (int i = 0; i < M; i ++)
    {
      v[i] = 0;
    }
    //v = null;
  }
}
