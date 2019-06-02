// Make it so there are collisions until an equilibrium is reached
//(nCollisions as the measure maybe?)

// then turn off collisions - maybe then turn on stats or have them on
// all the time to see distribution moving

// make the plots individual arrays so I can highlight each with object functions
// Look into matrix for partitioning
// find maxSpeedSq from the speed dist.
// make smaller plots for collisions and framerate
// Find density of particles in different regions of the box
// Have sum over all momentum to show it is 0?

int N = 10000000;
int averagingTime = 1; // how often to calculate statistics over draw loops
float dRad = 1;//0.03*depth; // R ~= depth/(4000*N)^(1/3) for 0.1% of total Volume to be particles
int M = 200;
int div = 200; // more bins leaves more dead room for collisions not to occour
int Nbins = div*div*div;

// the second dimension needs to have longer arrays around the expected velocity 
int binSecDim = floor(5*N/(Nbins));
// Taylor the length of the array to more efficiently fit the distribution?
int[][] bins = new int[Nbins][binSecDim]; // twice the expected number should be enough?
int[] binCount = new int [Nbins];

float[] xScale = new float[M];
float maxSp = 10;
float normFactor = 0.35*M/(N*pow(maxSp, float(1/3)));
float xStretch = 1;
int setFrameRate = 600;
float dDensity = 1;
float depth = 1000; // for 3D visuals

// Normalising factors to help plotting later
float maxPressure = 0;
float maxEnergy = 0;
int maxCollisions = 0;
float vTh;

// Maxwellian to match the plot to
float[] max = new float[M];

//float t = 0; // "time" as an x-axis
int count = 0; // counting variable to set how often we reset the statistical 0's

float piSpeedSq;
float sQspeedSum = 0;

int nCollisions = 0; // Number of collisions 

float energy;
float chi = 0;

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
int setupCount = 0;
int nEffuse = 0;

void setup()
{  
  size(1000, 1000);//, P3D);
  frameRate(setFrameRate);
  println("Step up");

  background(0);

  for (int i = 0; i < M; i ++)
  {
    v[i] = 0;
  }

  // setting scale for x axis
  for (int i = 0; i < M; i ++)
  {
    xScale[i] = xStretch*i*width/M;
  }

  // Make N particles with random positions and velocities
  // Add them to the ArrayList p
  for (int i = 0; i < N; i = i + 1) 
  {

    if (i % (N/100) == 0)
    {
      println("Setup is " + setupCount + "% done.");
      setupCount ++;
    }
    // Random Starting positions
    x = random(width);
    y = random(height);
    z = random(depth);

    vx = (i-N/2)*maxSp/N;

    p.add(new Particle(x, y, z, vx, vx, vx, dRad, dDensity));
  }
  println("Finished creating particles");

  updateDistribution();
  plot(v, 1, 255, 255, 255);
  maxWellian();
  plotMax();
  println("Finished Setup");
}

void draw()
{ 
  // These five functions are for 3D plotting of the actual gas 
  //background(0);
  //lights();
  //rotation();
  //borders();
  //axis();


  // reseting the colour of the particles after they have collided for 3D plot
  //for (int i = 0; i < N; i ++)
  //{
  //  p.get(i).col = color(255, 255, 255);
  //}

  // set running total of force to 0 at the start of each loop over averagingTime
  if (count == 0)
  {
    force = 0;
    //nCollisions = 0;
    sQspeedSum = 0;
  }
  count += 1;

  // Using the spatial bins to optimise the algorithm
  partition();

  // Run the collision detection for each of the bins
  for (int k = 0; k < Nbins; k ++)
  {

    // Looping through the ArrayList
    for (int i = 0; i < binCount[k] && i < binSecDim; i ++) 
    {

      // looping through the remaining particles for each particle to collide with
      for (int j = i + 1; j < binCount[k] && j < binSecDim; j ++) 
      {
        if (i != j)
        {
          // shortening the object name for the use in the loop
          pi = p.get(bins[k][i]);
          pj = p.get(bins[k][j]);

          // relative vector and scalar values for the displacements
          posDiff = pi.pos.copy().sub(pj.pos.copy());
          posDiffSq = magSq(posDiff);

          // Check to see if particles are close enough to interact mechanically
          if (posDiffSq <= 4*dRad*dRad)//(pi.rad + pj.rad)*(pi.rad + pj.rad))
          {
            collide(pi, pj);
            nCollisions ++;
          }
          //else
          //{
          // Apply any other forces

          //gravity(pi, pj);
          //Electric(pi, pj);
          //}
        }
      }
    }
  }

  // Different loop to magage particles to not have conflicts while doing calculations
  for (int i = 0; i < N; i = i + 1) 
  {
    Particle particle = p.get(i);      
    particle.edges();  // Container walls collisions
    particle.update();  // apply all forces and changes
    //p.get(i).show();  // Show the particles as spheres
  }

  //// Display statistics every averagingTime draw loops
  if (count >= averagingTime)
  {
    //background(0);
    count = 0;

    updateMaxVals();

    updateDistribution();
    calcThermSp();
    chiSq();
    displayParameters();
    stateVariablePlot();

    rainbow();
    plot(v, 0.3, 255 - colour, 0, colour);
    plotMax();

    nEffuse = 0;

    save("Maxwellian.jpg");
  }
}
