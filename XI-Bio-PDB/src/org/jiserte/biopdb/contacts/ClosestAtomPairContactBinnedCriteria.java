package org.jiserte.biopdb.contacts;

import pair.Pair;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jiserte.biopdb.structures.Atom;
import org.jiserte.biopdb.structures.Chain;
import org.jiserte.biopdb.structures.Residue;
import org.jiserte.biopdb.structures.SpacePoint;

/**
 * Is a criteria to decide if two residues are in a protein are in contact.
 * If the distance of the two closest atoms if two residues is less than a given
 * cutoff value then is assumed as contact.
 * 
 * @author javier
 *
 */
public class ClosestAtomPairContactBinnedCriteria extends ContactCriteria {

	////////////////////////////////////////////////////////////////////////////
	// Instance Variables
	private double cutoffDistance;
  private Map<Character, Chain> pdb;
  private boolean evaluated = false;
  private Map<Pair<Residue, Residue>, Pair<SpacePoint, SpacePoint>>
    evaluationResults;
	////////////////////////////////////////////////////////////////////////////

	
	////////////////////////////////////////////////////////////////////////////
	// Constructor
	/**
	 * Creates a new ClosestAtomPairContactCriteria with a given cut-off value.
	 * The distance is measured in Angstroms.
	 * @param cutoffDistance
	 */
	public ClosestAtomPairContactBinnedCriteria(double cutoffDistance) {
		this.setCutoffDistance(cutoffDistance);
	}
	////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////
	// Public Interface
	/**
	 * Check if two residues are in contact. If so, the closest atom pair is
	 * returned. Otherwise null is returned.
	 */
	@Override
	public Pair<SpacePoint, SpacePoint> areInContact(Residue currentFirstPoint,
			Residue currentSecondPoint) {

	  if (!this.evaluated) {
	    this.evaluate();
	  }
	  
	  
	  
	  Pair<Residue, Residue> testPair = getSortedResiduePair(currentFirstPoint,
	      currentSecondPoint);
	  
	  if (this.evaluationResults.containsKey(testPair)) {
	    return this.evaluationResults.get(testPair);
	  } else {
	    return null;
	  }
		
	}
	
	/**
	 * Search all residue pairs that are in contact
	 */
	private void evaluate() {
    Map<Integer,List<SpacePoint>> bins = new HashMap<>();
    Map<Pair<Residue,Residue>,Pair<SpacePoint,SpacePoint>> results = new HashMap<>();
    ////////////////////////////////////////////////////////////////////////////
    // Distribute atoms in bins according to the distance to the origin.  
    SpacePoint origin = new Atom();
    origin.setX(0);
    origin.setY(0);
    origin.setZ(0);
        for (Chain chain : pdb.values()) {
      for (Residue residue : chain.getResiduesCollection()) {
        for (SpacePoint point : residue.getSpacePointCollection()) {
          Integer binId = (int) (point.distanceTo(origin) / this.getCutoffDistance());
          if (!bins.containsKey(binId)) {
            bins.put(binId, new ArrayList<>());
          }
          bins.get(binId).add(point);
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////
        
    ////////////////////////////////////////////////////////////////////////////
    // Check all points in a bin
    ArrayList<Integer> orderedBins = new ArrayList<>();
    orderedBins.addAll(bins.keySet());
    Collections.sort(orderedBins);
    
    int min = orderedBins.get(0);
    int max = orderedBins.get(orderedBins.size()-1);
        
    
    for ( int binIndex = min+1; binIndex <= max ; binIndex++ ) {

      //////////////////////////////////////////////////////////////////////////
      // Creates a single list with the points of two adjacent bins.  
      List<SpacePoint> pointsInBin = new ArrayList<>();
      for (int i=0 ; i<2; i++) {
        if (bins.containsKey(binIndex-i)) {
          List<SpacePoint> morePoints = bins.get(binIndex-i);
          pointsInBin.addAll(morePoints);
        }
      }
      //////////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////////
      // Evaluate all pairs of atoms in two adjacent bins
      for (int index = 0; index < pointsInBin.size()-1; index++) {
        for (int index2 = index + 1; index2 < pointsInBin.size(); index2++) {
          
          
          SpacePoint atom1 = pointsInBin.get(index);
          SpacePoint atom2 = pointsInBin.get(index2);
          
          double distance = atom1.distanceTo(atom2);
          if ( distance <= this.getCutoffDistance()) {
            
            char chain1 = atom1.getChainIdentifier();
            char chain2 = atom2.getChainIdentifier();
            int resNumber1 = atom1.getResidueSequenceNumber(); 
            int resNumber2 = atom2.getResidueSequenceNumber();
            
            Residue residue1 = this.pdb.get(chain1).getResidues().get(resNumber1);
            Residue residue2 = this.pdb.get(chain2).getResidues().get(resNumber2);
            
            Pair<Residue,Residue> inContact = getSortedResiduePair(residue1,residue2 );
           
            if (results.containsKey(inContact)) {
              Pair<SpacePoint, SpacePoint> old = results.get(inContact);
              double oldDistance = old.getFirst().distanceTo(old.getSecond());
              
              if (distance<oldDistance) {
                results.put(inContact, new Pair<SpacePoint, SpacePoint>( 
                    atom1,  atom2));
              }
              
            } else {
              results.put(inContact, new Pair<SpacePoint, SpacePoint>( 
                  atom1,  atom2));
            }
          }
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////
    this.evaluated = true;
    this.evaluationResults = results;
    
  }
  /**
	 * Returns a value that indicates the cost of the operation
	 */
  @Override
  public int cost() {
    return 1000;
  }
	
	/**
	 * Gets the cut-off value to decide if two residues are in contact.
	 * @return
	 */
	public double getCutoffDistance() {
		return cutoffDistance;
	}
	/**
	 * Sets the cut-off value to decide if two residues are in contact.
	 * @return
	 */	
	public void setCutoffDistance(double cutoffDistance) {
		this.cutoffDistance = cutoffDistance;
	}
	

  @Override
  public boolean useDistance() {
    return true;
  }

  @Override
  public double getUsedDistance() {
    return this.getCutoffDistance();
  }
	////////////////////////////////////////////////////////////////////////////

  @Override
  public boolean requiresEntirePDB() {
    return true;
  }

  @Override
  public void setEntirePDB(Map<Character, Chain> pdb) {
    this.pdb = pdb;
  }

  @Override
  public boolean canProvideCandidates() {
    return true;
  }

  @Override
  public Set<Pair<Residue, Residue>> getCandidates() {
    if (! this.evaluated) {
      this.evaluate();
    }
    return this.evaluationResults.keySet();
    
  }
  
  private Pair<Residue, Residue> getSortedResiduePair(Residue res1, Residue res2) {
    
    char chainId1 = res1.getChain();
    char chainId2 = res2.getChain();
    
    if ( chainId1 < chainId2 ) {
      return new Pair<Residue,Residue>(res1,res2);
    } else 
    if (chainId2 < chainId1) {
      return new Pair<Residue,Residue>(res2,res1);
    } else {
      Integer resNum1 = res1.getResNumber();
      Integer resNum2 = res2.getResNumber();
      
      if ( resNum1 < resNum2 ) {
        return new Pair<Residue,Residue>(res1,res2);
      } else {
        return new Pair<Residue,Residue>(res2,res1);
      }
      
    }
  }

}
