class Particle
{
  PVector pos, vel, acc;
  float rad, density, mass, charge, vol, col;

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

    col = map(density, dDensity, 5 * dDensity, 255, 0);
  }

  void update()
  {
    vel.add(acc);
    pos.add(vel);
    acc.set(zero.copy());
    //vel.limit(maxVel);
  }

  void show()
  {
    fill(col, col, 255);
    noStroke();
    pushMatrix();
    translate(pos.x, pos.y, pos.z);
    sphere(rad);
    popMatrix();
  }

  void addForce(PVector force)
  {
    acc.add(force.copy().div(mass));
  }

  void edges()
  {
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
  return num * num * num;
}


//void axis()
//{
//  float lineLen = 100;
//  strokeWeight(4);
//  translate(width/2, height/2);

//  stroke(255, 0, 0); // RED X axis
//  line(0, lineLen, 0, 0, 0, 0);

//  stroke(0, 255, 0); // GREEN Y axis
//  line(0, 0, 0, lineLen, 0, 0);

//  stroke(0, 0, 255); // BLUE Z axis
//  line(0, 0, 0, 0, 0, lineLen);

//  translate(-width/2, -height/2);
//}



void rotation() // fix
{
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
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

void borders()
{
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
  for (int i = 0; i < div*div; i ++)
  {
    for (int j = 0; j < binSecDim; j ++)
    {
      bins[i][j] = 0;
    }
  }
}

void partition()
{
  resetBins();

  for (int i = 0; i < Nbins; i ++)
  {
    binCount[i] = 0;
  }

  for (int j = 0; j < N; j ++)
  {
    int x = floor((p.get(j).pos.x / width) * div);
    int y = floor((p.get(j).pos.y / height) * div);

    int binIndex = x + y * div;// - div - 1;

    //println("x: " + x + ", y: " + y);
    if (binIndex >= 0 && binIndex < Nbins)
    {
      binCount[binIndex] += 1;
      if (binCount[binIndex] < binSecDim)
      {
        bins[binIndex][binCount[binIndex]] = j;
      }
    }
    // put the index of the particle in the bin that it is in
    //
  }
}

void plot(int[] array, float weight, int rC, int bC, int gC)
{
  stroke(rC, bC, gC);
  strokeWeight(weight);

  for (int i = 0; i < array.length - 1; i ++)
  {
    //point(1*i*width/M, height*(1 - normFactor*v[i]));

    line(5 + xScale[i], height*(1 - normFactor*array[i]), 
      5+xScale[i + 1], height*(1 - normFactor*array[i + 1]));
  }
}

void updateDistribution()
{
  for (int i = 0; i < N; i = i + 1)
  {
    //p.get(i).show();
    piSpeedSq = magSq(p.get(i).vel);
    //println(piSpeedSq);
    //sQspeedSum += piSpeedSq;

    index = floor(0.1*M * piSpeedSq/maxSpSq);

    if (index < M)
    {
      v[index] += 1;
    }
  }
}

void stateVariablePlot()
{
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
  println("Collisions: " + float(nCollisions)/averagingTime);
  //println("Average Speed Squared: " + energy);
  //println("N: " + N);
  //println();
}

void updateMaxVals()
{
  //    if (pressure > maxPressure)
  //    { maxPressure = pressure;}
  //    if (energy > maxEnergy)
  //    {maxEnergy = energy;}
  if (nCollisions > maxCollisions)
  {
    maxCollisions = nCollisions;
  }
}
