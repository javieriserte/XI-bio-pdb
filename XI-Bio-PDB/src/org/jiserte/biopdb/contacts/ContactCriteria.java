package org.jiserte.biopdb.contacts;

import pair.Pair;
import org.jiserte.biopdb.structures.Residue;
import org.jiserte.biopdb.structures.SpacePoint;

public abstract class ContactCriteria {

  public abstract Pair<SpacePoint, SpacePoint> areInContact(
      Residue currentFirstPoint, Residue currentSecondPoint);
  
  public abstract boolean useDistance();
  
  public abstract double getUsedDistance();

  public abstract int cost();

}
