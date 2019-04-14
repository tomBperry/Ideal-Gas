// Make it so there are collisions until an equilibrium is reached
//(nCollisions as the measure maybe?)

// then turn off collisions - maybe then turn on stats or have them on
// all the time to see distribution moving

// make the plot root(v) on x-axis
// increase scale on y-axis
// decrease the colour change per draw loop
// make the plots individual arrays so I can highlight each with object functions
// Look into matrix for partitioning


int N = 500000;
int averagingTime = 1; // how often to calculate statistics over draw loops
float dRad = 3; // R ~= (V/(4000*N))^(1/3) for 0.1% of total Volume to be particles
int M = 200;

// Must be a squared factor of N
int div = 40; // more bins leaves more dead room for collisions not to occour

// the second dimension needs to have longer arrays around the expected velocity 
int[][] bins = new int[div*div][2*N/div*div]; // twice the expected number should be enough?
int[] binCount = new int [div*div];


float normFactor = 0.3*float(M)/float(N);
int setFrameRate = 1000;
float dDensity = 1;
int randVelAdd = 1; // max random initial velocity in each direction (x, y, z)
//float maxVel = 2;
float depth = 1000; // for 3D visuals

float maxSpSq = 4;

// Normalising factors to help plotting later
float maxPressure = 0;
float maxEnergy = 0;
int maxCollisions = 0;

float t = 0; // "time" as an x-axis
int count = 0; // counting variable to set how often we reset the statistical 0's

float piSpeedSq;
float sQspeedSum = 0;

int nCollisions = 0; // Number of collisions 

float energy;

int colour = 0;


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
//int  
PVector posDiff;

int index;

void setup()
{  
  size(1000, 1000);//, P3D);
  frameRate(setFrameRate);

  background(0);

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
    energy += vx*vx + vy*vy + vz*vz;
    //    size = p.size();
  }
  energy /= N;

  for (int i = 0; i < N; i = i + 1)
  {
    //p.get(i).show();
    piSpeedSq = magSq(p.get(i).vel);
    sQspeedSum += piSpeedSq;

    index = floor(M * piSpeedSq / maxSpSq);
    //println(index);
    //if (index < M)
    //{
    v[index] += 1;
    //}
  }

  for (int i = 0; i < M-1; i ++)
  {
    //if (v[i] != 0)
    //{
    //println("count " + i + ": " + v[i]);
    //}
    //point(1*i*width/M, height*(1 - normFactor*v[i]/averagingTime));
    stroke(255);
    strokeWeight(0.5);
    line(1*i*width/M, height*(1 - normFactor*v[i]/averagingTime), 
      1*(i+1)*width/M, height*(1 - normFactor*v[i + 1]/averagingTime));
    // make this a line graph
  }
  println("Hello");
}



void draw()
{ 
  //background(0);
  //lights();


  //rotation();
  //borders();


  // set running total of force to 0 at the start of each loop

  if (count == 0)
  {
    force = 0;
    nCollisions = 0;
    sQspeedSum = 0;
  }
  count += 1;
  //println(count);

  // Loop is how many times the system evolves per draw loop
  //for (int loop = 1; loop > 0; loop = loop - 1)
  //{
  partition();

  for (int k = 0; k < div*div; k ++)
  {

    // Looping through the ArrayList
    for (int i = binCount[k]; i >= 0; i = i - 1) 
    {

      // looping through the remaining particles for each particle
      for (int j = i-1; j >= 0; j = j - 1) 
      {
        if (i != j)
        {
          // shortening the object name for the use in the loop
          pi = p.get(bins[k][i]);
          pj = p.get(bins[k][j]);

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

  for (int i = 0; i < M; i ++)
  {
    v[i] = 0;
  }

  for (int i = 0; i < N; i = i + 1)
  {
    //p.get(i).show();
    piSpeedSq = magSq(p.get(i).vel);
    sQspeedSum += piSpeedSq;

    index = floor(M * piSpeedSq / maxSpSq);
    //println(index);
    if (index < M)
    {
      v[index] += 1;
    }
  }


  // Display statistics every time that averagingTime loops
  if (count >= averagingTime)
  {
    //background(0);
    count = 0;



    // Show all the particles in their new positons 
    //for (int i = 0; i < N; i = i + 1)
    //{
    //  //p.get(i).show();
    //  piSpeedSq = magSq(p.get(i).vel);
    //  sQspeedSum += piSpeedSq;

    //  int index = floor(M * piSpeedSq / maxSpSq);
    //  //println(index);
    //  v[index] += 1;
    //}

    if (pressure > maxPressure)
    {
      maxPressure = pressure;
    }
    if (energy > maxEnergy)
    {
      maxEnergy = energy;
    }
    if (nCollisions > maxCollisions)
    {
      maxCollisions = nCollisions;
    }

    // Calculating the pressure from the wall collisions
    // Could save this data and plot it - Thermodynamic limit
    pressure = force;//(6*width*height)*1000;
    energy = sQspeedSum;



    //stroke(255, 0, 0); // Red Pressure
    //point(t, height*(1-pressure/maxPressure));

    //stroke(0, 255, 0); // Green Energy
    //point(t, height*(1-energy/maxEnergy));

    //stroke(0, 0, 255); // Blue Collisions
    //point(t, height*(1-float(nCollisions)/float(maxCollisions)));

    println("FrameRate: " + frameRate);
    //println("Pressure: " + pressure);
    println("Collisions: " + float(nCollisions)/averagingTime);
    println("Average Speed Squared: " + energy);
    //println("N: " + N);
    //println();

    if (colour >= 255)
    {
      colour %= 255;
    } else
    {
      colour += 10;
    }
    stroke(255-colour, 0, colour);
    strokeWeight(0.5);


    for (int i = 0; i < M - 1; i ++)
    {
      //if (v[i] != 0)
      //{
      //println("count " + i + ": " + v[i]);
      //}
      //point(1*i*width/M, height*(1 - normFactor*v[i]/averagingTime));
      line(1*i*width/M, height*(1 - normFactor*v[i]/averagingTime), 
        1*(i+1)*width/M, height*(1 - normFactor*v[i + 1]/averagingTime));
      // make this a line graph
    }

    // Add a best fit function

    save("Maxwellian.jpg");
    //could save the images formed?? Save all in a sequence in a folder?


    if (t > width)
    {
      t = 0;
    }
    t += 1;

    //v = null;
  }
  //}
}
