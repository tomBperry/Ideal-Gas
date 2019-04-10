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
    vel.limit(maxVel);
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


void axis()
{
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
