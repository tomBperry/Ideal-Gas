// Make it so there are collisions until an equilibrium is reached
//(nCollisions as the measure maybe?)

// then turn off collisions - maybe then turn on stats or have them on
// all the time to see distribution moving

// make the plot root(v) on x-axis
// increase scale on y-axis
// decrease the colour change per draw loop
// make the plots individual arrays so I can highlight each with object functions
// Look into matrix for partitioning
// find maxSpeedSq from the speed dist.
// make smaller plots for collisions and framerate 

int N = 100000;
int averagingTime = 1; // how often to calculate statistics over draw loops

float depth = 1000; // for 3D visuals
float dRad = 2;//0.03*depth; // R ~= depth/(4000*N)^(1/3) for 0.1% of total Volume to be particles
int M = 1000;

// Must be a squared factor of N
int div = 5; // more bins leaves more dead room for collisions not to occour
int Nbins = div*div;
// the second dimension needs to have longer arrays around the expected velocity 
int binSecDim = floor(2*N/(Nbins));
int[][] bins = new int[Nbins][binSecDim]; // twice the expected number should be enough?
int[] binCount = new int [Nbins];

float[] xScale = new float[M];


float normFactor = 0.001;//0.float(M)/(float(N)*averagingTime);
float xStretch = 30;
int setFrameRate = 1000;
float dDensity = 1;
int randVelAdd = 1; // max random initial velocity in each direction (x, y, z)
//float maxVel = 2;


float maxSpSq = 4;

// Normalising factors to help plotting later
float maxPressure = 0;
float maxEnergy = 0;
int maxCollisions = 0;

//float t = 0; // "time" as an x-axis
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

  // setting scale for x axis
  for (int i = 0; i < M; i ++)
  {
    xScale[i] = xStretch*pow(i*width/M, 0.5);
  }

  // Make N particles with random positions and velocities
  // Add them to the ArrayList p
  for (int i = 0; i < N; i = i + 1) 
  {
    x = random(width);
    y = random(height);
    z = random(depth);

    vx = pow(i*maxSpSq/N, 0.5);

    p.add(new Particle(x, y, z, vx, vx, vx, dRad, dDensity));
    //energy += vx*vx + vy*vy + vz*vz;
    //    size = p.size();
  }
  //energy /= N;

  updateDistribution();
  plot(v, 1, 255, 255, 255);
}

void draw()
{ 
  //background(0);
  //lights();
  //rotation();
  //borders();


  // set running total of force to 0 at the start of each loop over averagingTime
  if (count == 0)
  {
    //force = 0;
    nCollisions = 0;
    //sQspeedSum = 0;
  }
  count += 1;

  partition();

  for (int k = 0; k < Nbins; k ++)
  {

    // Looping through the ArrayList
    for (int i = binCount[k]; i >= 0 && i < binSecDim; i = i - 1) 
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
          if (posDiffSq <= 12*dRad*dRad)//(pi.rad + pj.rad)*(pi.rad + pj.rad))
          {
            collide(pi, pj);
            nCollisions ++;
          } else
          {
            // Apply any other forces

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
      particle.edges();  // Container walls collisions
      particle.update();  // apply all forces and changes
    }
  }

  for (int i = 0; i < M; i ++)
  {
    v[i] = 0;
  }

  updateDistribution();

  // Display statistics every averagingTime draw loops
  if (count >= averagingTime)
  {
    //background(0);
    count = 0;

    updateMaxVals();

    stateVariablePlot();
    plot(v, 0.1, 255 - colour, 0, colour); 

    //Calculating the pressure from the wall collisions
    // Could save this data and plot it - Thermodynamic limit
    //pressure = force;//(6*width*height)*1000;
    //energy = sQspeedSum;

    save("Maxwellian.jpg");
    //could save the images formed?? Save all in a sequence in a folder?
  }
}
