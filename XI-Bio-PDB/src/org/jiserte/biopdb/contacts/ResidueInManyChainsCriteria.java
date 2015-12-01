package org.jiserte.biopdb.contacts;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.jiserte.biopdb.structures.Chain;
import org.jiserte.biopdb.structures.Residue;

import pair.Pair;

/**
 * This criteria allows that two residue belongs to a set of predefined chains
 * 
 * @author javier
 *
 */
public class ResidueInManyChainsCriteria extends ChainCriteria {
  // ///////////////////////////////////////////////////////////////////////////
  // Instance Variables
  private Set<Character> allowedChains;

  // /////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////
  // Constructor
  public ResidueInManyChainsCriteria(Set<Character> allowedChains) {
    this.setAllowedChains(allowedChains);
  }

  public ResidueInManyChainsCriteria(String allowedChains) {
    Set<Character> allowed = new HashSet<>();
    for (Character c : allowedChains.toUpperCase().toCharArray()) {
      allowed.add(c);
    }
    this.setAllowedChains(allowed);
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////
  // Public Interface
  /**
   * @return the allowedChains
   */
  public Set<Character> getAllowedChains() {
    return allowedChains;
  }

  /**
   * @param allowedChains
   *          the allowedChains to set
   */
  public void setAllowedChains(Set<Character> allowedChains) {
    this.allowedChains = allowedChains;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////
  // Protected methods
  @Override
  protected boolean inContact(Residue currentFirstPoint,
      Residue currentSecondPoint) {
    return this.getAllowedChains().contains(currentFirstPoint.getChain()) && 
           this.getAllowedChains().contains(currentSecondPoint.getChain());
  }
  // ///////////////////////////////////////////////////////////////////////////

  @Override
  public boolean useDistance() {
    return false;
  }

  @Override
  public double getUsedDistance() {
    return 0;
  }

  @Override
  public boolean requiresEntirePDB() {
    return false;
  }

  @Override
  public void setEntirePDB(Map<Character, Chain> pdb) {}

  @Override
  public boolean canProvideCandidates() {
    return false;
  }

  @Override
  public Set<Pair<Residue, Residue>> getCandidates() {
    return new HashSet<>();
  }
}
