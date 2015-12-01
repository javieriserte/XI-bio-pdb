package org.jiserte.biopdb.contacts;

import pair.Pair;

import java.util.Map;
import java.util.Set;

import org.jiserte.biopdb.structures.Chain;
import org.jiserte.biopdb.structures.Residue;
import org.jiserte.biopdb.structures.SpacePoint;

public abstract class ContactCriteria {

  public abstract Pair<SpacePoint, SpacePoint> areInContact(
      Residue currentFirstPoint, Residue currentSecondPoint);
  
  public abstract boolean useDistance();
  
  public abstract double getUsedDistance();

  public abstract int cost();
  
  public abstract boolean requiresEntirePDB();
  
  public abstract void setEntirePDB(Map<Character, Chain> pdb);
  
  public abstract boolean canProvideCandidates();
  
  public abstract Set<Pair<Residue,Residue>> getCandidates();

}
