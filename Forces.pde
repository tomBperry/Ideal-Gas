void gravity(Particle particle1, Particle particle2)
{
  // Graviational force between two particles
  // Inputs: two particles (1 and 2)
  float G = 0.1; // 0.001 for sun and two particles  

  float gravMag = (G * particle1.mass * particle2.mass)/posDiffSq;
  PVector gravForce = posDiff.copy().setMag(gravMag);

  particle2.addForce(gravForce.copy());
  particle1.addForce(gravForce.copy().mult(-1));
}


void Electric(Particle particle1, Particle particle2)
{
  // Electirc force between two particles
  // Inputs: two particles (1 and 2)
  float Ec = -0.1; // 0.001 for sun and two particles  

  float elMag = (Ec * particle1.charge * particle2.charge)/posDiffSq;
  PVector elForce = posDiff.copy().setMag(elMag);

  particle2.addForce(elForce.copy());
  particle1.addForce(elForce.copy().mult(-1));
}

void collide(Particle p1, Particle p2)
{
  // Elastic binary collisons between two particles
  //Inputs: two particles (1 and 2)
  PVector v1v2Diff = p1.vel.copy().sub(p2.vel.copy());
  float vDiff = v1v2Diff.mag();

  PVector vVcm = p1.vel.copy();
  vVcm.add(p2.vel.copy());
  vVcm.div(2);

  PVector v1Temp = posDiff.copy().setMag(vDiff * 0.5);
  PVector v2Temp = posDiff.copy().setMag(-vDiff * 0.5);

  p1.vel.set(v1Temp.add(vVcm));
  p2.vel.set(v2Temp.add(vVcm));

  // Change their colour to red upon collision
  //p1.col = color(255, 0, 0);
  //p2.col = color(255, 0, 0);
}
