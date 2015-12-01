package org.jiserte.biopdb.contacts;

import annotations.NeverUsed;
import pair.Pair;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.jiserte.biopdb.structures.Chain;
import org.jiserte.biopdb.structures.Residue;


/**
 * This class represents a criteria to decide if two residues are in contact. In
 * this criteria, they must be in a given chain.
 * 
 * @author javier iserte
 *
 */
@NeverUsed
public class ResidueInChainCriteria extends ChainCriteria {
  // ///////////////////////////////////////////////////////////////////////////
  // Instance variables
  private char chain;
  // ///////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////
  // Constructor
  public ResidueInChainCriteria(char chain) {
    this.setChain(chain);
  }
  public ResidueInChainCriteria() {
    this.setChain('A');
  }
  // ///////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////  
  // Public interface
  /**
   * @return the chain
   */
  public char getChain() {
    return chain;
  }

  /**
   * @param chain
   *          the chain to set
   */
  public void setChain(char chain) {
    this.chain = chain;
  }
  // ///////////////////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////////////////
  // Protected interface
  @Override
  protected boolean inContact(Residue currentFirstPoint,
      Residue currentSecondPoint) {
    return currentFirstPoint.getChain() == this.getChain()
        && currentSecondPoint.getChain() == this.getChain();
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
