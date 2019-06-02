class Particle
{
  PVector pos, vel, acc;
  float rad, density, mass, charge, vol;
  color col;

  Particle(float x0, float y0, float z0, float vx0, float vy0, float vz0, float radius, float density)      
  {
    pos = new PVector(x0, y0, z0);
    vel = new PVector(vx0, vy0, vz0);
    acc = new PVector(0, 0, 0);

    rad = radius;
    vol = (4/3) * PI * cube(rad);
    this.density = density; 
    mass = vol * density;

    charge = 1;

    col = color(255, 255, 255);
  }

  void update()
  {
    // Apply the kinematics to the particles
    vel.add(acc);
    pos.add(vel);
    acc.set(zero.copy());
    //vel.limit(maxVel);
  }

  void show()
  {
    // Draw the particles as spheres at their positions in their colour
    fill(col);
    noStroke();
    pushMatrix();
    translate(pos.x, pos.y, pos.z);
    sphere(rad);
    popMatrix();
  }

  void addForce(PVector force)
  {
    // Add force to the acceleration through Newton 2
    acc.add(force.copy().div(mass));
  }

  void edges()
  {
    // If a particle hits a wall then reverse its velocity in that direction
    // if the particle is too far this may cause it to freeze in place

    // Find the force taken to do this and add it to the total

    if (this.pos.x > width)
    {
      this.vel.x *= -1;
      force -= 2*this.mass*this.vel.x;
    }

    if (this.pos.y > height)
    {
      this.vel.y *= -1;
      force -= 2*this.mass*this.vel.y;
    }

    if (this.pos.z > depth)
    {
      this.vel.z *= -1;
      force -= 2*this.mass*this.vel.z;
    }

    if (this.pos.x < 0)
    {
      this.vel.x *= -1;
      force += 2*this.mass*this.vel.x;
    } 

    if (this.pos.y < 0)
    {
      this.vel.y *= -1;
      force += 2*this.mass*this.vel.y;
    }

    if (this.pos.z < 0)
    {
      this.vel.z *= -1;
      force += 2*this.mass*this.vel.z;
    }
  }
}


float cube(float num)
{
  // Simple function to cube a float
  return num * num * num;
}


void axis()
{
  // Draw a caresian reference axis

  float lineLen = 100;
  strokeWeight(4);
  translate(width/2, height/2);

  stroke(255, 0, 0); // RED X axis
  line(0, lineLen, 0, 0, 0, 0);

  stroke(0, 255, 0); // GREEN Y axis
  line(0, 0, 0, lineLen, 0, 0);

  stroke(0, 0, 255); // BLUE Z axis
  line(0, 0, 0, 0, 0, lineLen);

  translate(-width/2, -height/2);
}



void rotation() // fix
{
  // rotate the 3D space for better viewing

  if (mousePressed)
  {
    rotAngleX = map(mouseY, 0, height, 0, TWO_PI)%TWO_PI;
    rotAngleY = map(mouseX, 0, width, 0, TWO_PI)%TWO_PI;
  }
  //println("rotAngleX: " + rotAngleX);
  //println("rotAngleY: " + rotAngleY);


  float x = width/2;
  float y = height/2;
  float z = depth/2;


  translate(x, y, 0);
  //translate(0, 0, z);
  rotateX(rotAngleX);
  //translate(0, -y, -z);

  //translate(x, 0, z);
  rotateY(rotAngleY);
  scale(0.4);
}


float magSq(PVector v)
{
  // Find the magnitude squared of a vector
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

void borders()
{
  // For 3D plotting draw a line box to enclose the particles
  strokeWeight(3);
  stroke(255);
  line(0, 0, 0, width, 0, 0);
  line(0, 0, 0, 0, height, 0);
  line(0, 0, 0, 0, 0, depth);
  line(width, 0, 0, width, height, 0);
  line(width, 0, 0, width, 0, depth);
  line(0, 0, depth, width, 0, depth);  
  line(0, 0, depth, 0, height, depth);
  line(0, height, 0, width, height, 0);
  line(0, height, 0, 0, height, depth);
  line(0, height, depth, width, height, depth);
  line(width, height, 0, width, height, depth);
  line(width, 0, depth, width, height, depth);
}

void resetBins()
{
  // Reset the arrays to 0 for each loop
  for (int i = 0; i < Nbins; i ++)
  {
    for (int j = 0; j < binSecDim; j ++)
    {
      bins[i][j] = 0;
    }
  }

  for (int i = 0; i < Nbins; i ++)
  {
    binCount[i] = 0;
  }
}

void partition()
{
  // Separate the particles into spacial bins depending on their position
  resetBins();

  for (int j = 0; j < N; j ++)
  {
    // Find which bin the particle should be in
    int x = floor((p.get(j).pos.x / width) * div);
    int y = floor((p.get(j).pos.y / height) * div);
    int z = floor((p.get(j).pos.z / depth) * div);
    int binIndex = x + y * div + z * div*div; 

    // Effusions are for the moment representing how many particles are in each bin
    if (binIndex == 0)
    {
      nEffuse += 1;
    }

    //println("x: " + x + ", y: " + y);
    if (binIndex >= 0 && binIndex < Nbins)
    {
      // if that bin isn't full (max set for optimisation)
      if (binCount[binIndex] < binSecDim)
      {
        // specifiy the index to that bin
        bins[binIndex][binCount[binIndex]] = j;
      }

      // incriment the number of particles in that bin
      binCount[binIndex] += 1;
    }
  }
}

void plot(int[] array, float weight, int rC, int bC, int gC)
{
  // Plot the speed distribtution with certain colours and line thickness
  stroke(rC, bC, gC);
  strokeWeight(weight);

  for (int i = 0; i < array.length - 1; i ++)
  {
    //point(1*i*width/M, height*(1 - normFactor*v[i]));

    line(xScale[i], height*(1 - normFactor*array[i]), 
      xScale[i + 1], height*(1 - normFactor*array[i + 1]));
  }
}

void updateDistribution()
{
  // For each frame, sort the particles into M speed bins for plotting
   
  // Reseting the array to 0
  for (int i = 0; i < M; i ++)
  {
    v[i] = 0;
  }
  
  for (int i = 0; i < N; i = i + 1)
  {
    piSpeedSq = p.get(i).vel.magSq();
    sQspeedSum += piSpeedSq;

    index = floor(0.8*sqrt(piSpeedSq)*M/maxSp);

    if (index < M)
    {
      v[index] += 1;
    }
  }
}

void stateVariablePlot()
{
  // print to the console various statistical properties of the simulation
  
  //stroke(255, 0, 0); // Red Pressure
  //point(t, height*(1-pressure/maxPressure));
  //stroke(0, 255, 0); // Green Energy
  //point(t, height*(1-energy/maxEnergy));
  //stroke(0, 0, 255); // Blue Collisions
  //point(t, height*(1-float(nCollisions)/float(maxCollisions)));

  //if (t > width)
  //{
  //  t = 0;
  //}
  //t += 1;


  println("FrameRate: " + frameRate);
  //println("Pressure: " + pressure);
  //println("Collisions: " + float(nCollisions));
  println("Collisions/Particle: " + float(nCollisions)/float(N));
  println("Energy/N: " + 0.5*p.get(0).mass*sQspeedSum/N);
  println("Effusions: " + nEffuse/averagingTime);
  //println("N: " + N);
  //println();
}

void updateMaxVals()
{
  // These values are used as normalisation constants for different plots if needed
  
  //    if (pressure > maxPressure)
  //    { maxPressure = pressure;}
  //    if (energy > maxEnergy)
  //    {maxEnergy = energy;}
  if (nCollisions > maxCollisions)
  {
    maxCollisions = nCollisions;
  }
}

void calcThermSp()
{
  // Find the mode of the speed distribution and set this as the thermal speed
  
  int highestIndex = 0;
  for (int i = 0; i < M; i ++)
  {
    if (v[i] > highestIndex)
    {
      highestIndex = i;
    }
  }

  vTh = highestIndex*maxSp/M;
}

void displayParameters()
{
  // Display important statistics in the top right of the plot
  
  int xOffset = 300;
  int tSize = 20;
  fill(0);
  stroke(0);
  rect(width - xOffset, 0, width, tSize*4);

  stroke(255);
  fill(255);
  textSize(20);
  text("N: " + str(N), width - xOffset, tSize);
  text("M: " + str(M), width - xOffset, 2*tSize);
  text("Collisions/Particle: " + float(nCollisions)/float(N), width - xOffset, 3*tSize);
  text("Chi Squared: " + chi, width - xOffset, 4*tSize);
}

void rainbow()
{
  // incriment the colour of the speed distribution plot
  
  colour += 1;
  if (colour > 255)
  {
    colour = 0;
  }
}

void maxWellian()
{
  // Make an array of the funcional form of a Maxwellian distribution to compare the derived one against
  
  float k = 230;
  float maxFactor = 440;
  for (int i = 0; i < M; i++)
  {
    float x = xScale[i];
    max[i] = maxFactor*x*x*exp(-x*x/(2*k*k));
    max[i] /= (k*k*k*PI*0.5);
    println(max[i]);
  }
}

void plotMax()
{
  // Plot the functional form of the Maxwellian
  
  stroke(255);
  strokeWeight(1);

  for (int i = 0; i < M - 1; i ++)
  {
    //point(1*i*width/M, height*(1 - normFactor*v[i]));

    line(xScale[i], height*(1 - max[i]), 
      xScale[i + 1], height*(1 - max[i + 1]));
  }
}

void chiSq()
{
  // Statistical test for the fitness of the two plots  
  
  // Sum of difference between observed and expected value all squared over the expected value
  chi = 0;
  for (int i = 0; i < M - 1; i ++)
  {
    if (max[i] != 0)
    {
      chi += pow((max[i] - normFactor*v[i]), 2)/max[i];
    }
  }
  println("Chi Squared value: " + chi);
}
