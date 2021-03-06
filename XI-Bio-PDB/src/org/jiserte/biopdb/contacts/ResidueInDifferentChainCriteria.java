package org.jiserte.biopdb.contacts;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.jiserte.biopdb.structures.Chain;
import org.jiserte.biopdb.structures.Residue;

import pair.Pair;

/**
 * This class represents a criteria to decide if two residues are in contact. In
 * this criteria, they must be in different chains.
 * 
 * @author javier iserte
 *
 */
public class ResidueInDifferentChainCriteria extends ChainCriteria {
  // ///////////////////////////////////////////////////////////////////////////
  // Instance variables
  private Set<Character> acceptedChains;

  // ///////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////
  // Constructors
  public ResidueInDifferentChainCriteria(Set<Character> acceptedChains) {
    this.setAcceptedChains(acceptedChains);
  }

  public ResidueInDifferentChainCriteria() {
    this.setAcceptedChains(new HashSet<Character>());
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////
  // Public interface
  /**
   * @return the acceptedChains
   */
  public Set<Character> getAcceptedChains() {
    return acceptedChains;
  }

  /**
   * @param acceptedChains
   *          the acceptedChains to set
   */
  public void setAcceptedChains(Set<Character> acceptedChains) {
    this.acceptedChains = acceptedChains;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////
  // Protected methods
  @Override
  protected boolean inContact(Residue currentFirstPoint,
      Residue currentSecondPoint) {
    Character chainFirst = currentFirstPoint.getChain();
    Character chainSecond = currentSecondPoint.getChain();
    return this.acceptedChains.contains(chainFirst)
        && this.acceptedChains.contains(chainSecond)
        && chainFirst != chainSecond;
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
