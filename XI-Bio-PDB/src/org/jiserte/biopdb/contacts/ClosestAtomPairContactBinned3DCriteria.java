package org.jiserte.biopdb.contacts;

import pair.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
public class ClosestAtomPairContactBinned3DCriteria extends ContactCriteria {

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
	public ClosestAtomPairContactBinned3DCriteria(double cutoffDistance) {
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
	  
	  Pair<Residue, Residue> testPair = this.getSortedResiduePair(currentFirstPoint,
	      currentSecondPoint);
	  
	  if (testPair != null && this.evaluationResults.containsKey(testPair)) {
	    return this.evaluationResults.get(testPair);
	  } else {
	    return null;
	  }
		
	}
	
	/**
	 * Search all residue pairs that are in contact
	 */
	private void evaluate() {
    Map<List<Integer>,List<SpacePoint>> bins = new HashMap<>();
    Map<Pair<Residue,Residue>,Pair<SpacePoint,SpacePoint>> results = new HashMap<>();

    ////////////////////////////////////////////////////////////////////////////
    // Variables to store the borders of each dimension.  
    int minBinX=Integer.MAX_VALUE;
    int maxBinX=Integer.MIN_VALUE;
    int minBinY=Integer.MAX_VALUE;
    int maxBinY=Integer.MIN_VALUE;
    int minBinZ=Integer.MAX_VALUE;
    int maxBinZ=Integer.MIN_VALUE;
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Distribute atoms in bins according to the distance to the origin.  
    for (Chain chain : pdb.values()) {
      for (Residue residue : chain.getResiduesCollection()) {
        for (SpacePoint point : residue.getSpacePointCollection()) {
          
          Integer binIdX = (int) ( Math.ceil(point.getX()
              / this.getCutoffDistance()));
          Integer binIdY = (int) (Math.ceil(point.getY()
              / this.getCutoffDistance()));
          Integer binIdZ = (int) (Math.ceil(point.getZ()
              / this.getCutoffDistance()));
          
          minBinX=Math.min(minBinX, binIdX);
          maxBinX=Math.max(maxBinX, binIdX);
          minBinY=Math.min(minBinY, binIdY);
          maxBinY=Math.max(maxBinY, binIdY);
          minBinZ=Math.min(minBinZ, binIdZ);
          maxBinZ=Math.max(maxBinZ, binIdZ);
          
          List<Integer> binId = createTriplet(binIdX, binIdY, binIdZ);
          
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
    for (int i = minBinX+1 ; i <= maxBinX ; i++) {
      for (int j = minBinY+1 ; j <= maxBinY ; j++) {
        for (int k = minBinZ+1 ; k <= maxBinZ ; k++) {
          
          //////////////////////////////////////////////////////////////////////
          // Creates a single list with the points of ALL adjacent bins.  
          List<SpacePoint> pointsInBin = new ArrayList<>();
          for (int ii=0 ; ii<2; ii++) {
            for (int jj=0 ; jj<2; jj++) {
              for (int kk=0 ; kk<2; kk++) {
                
                List<Integer> index = createTriplet(i-ii, j-jj, k-kk); 
                    
                if (bins.containsKey(index)) {
                  List<SpacePoint> morePoints = bins.get(index);
                  pointsInBin.addAll(morePoints);
                }
              }
            }
          }
          //////////////////////////////////////////////////////////////////////
          
          //////////////////////////////////////////////////////////////////////
          // Evaluate all pairs of atoms in All adjacent bins
          this.evaluteAtomsInBin(results, pointsInBin);
          //////////////////////////////////////////////////////////////////////
        }
      }
    }
    ////////////////////////////////////////////////////////////////////////////
    this.evaluated = true;
    this.evaluationResults = results;
    
  }

  private void evaluteAtomsInBin(
      Map<Pair<Residue, Residue>, Pair<SpacePoint, SpacePoint>> results,
      List<SpacePoint> pointsInBin) {
    for (int index = 0; index < pointsInBin.size()-1; index++) {
      for (int index2 = index + 1; index2 < pointsInBin.size(); 
          index2++) {
        
        SpacePoint atom1 = pointsInBin.get(index);
        SpacePoint atom2 = pointsInBin.get(index2);
        
        double distance = atom1.distanceTo(atom2);
        if ( distance <= this.getCutoffDistance()) {
          
          char chain1 = atom1.getChainIdentifier();
          char chain2 = atom2.getChainIdentifier();
          int resNumber1 = atom1.getResidueSequenceNumber(); 
          int resNumber2 = atom2.getResidueSequenceNumber();
          
          Residue residue1 = this.pdb.get(chain1).getResidues().get(
              resNumber1);
          Residue residue2 = this.pdb.get(chain2).getResidues().get(
              resNumber2);
          
          Pair<Residue,Residue> inContact = this.getSortedResiduePair(
              residue1,residue2 );
          if (inContact != null) {           
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
  }

  private List<Integer> createTriplet(Integer binIdX, Integer binIdY,
      Integer binIdZ) {
    List<Integer> binId = new ArrayList<>();
    binId.add(binIdX);
    binId.add(binIdY);
    binId.add(binIdZ);
    return binId;
  }
  /**
	 * Returns a value that indicates the cost of the operation
	 */
  @Override
  public int cost() {
    return 10;
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
    Set<Pair<Residue, Residue>> keySet = this.evaluationResults.keySet();
    return keySet;
    
  }
  
  private Pair<Residue, Residue> getSortedResiduePair(Residue res1, Residue res2) {
    
    char chainId1 = res1.getChain();
    char chainId2 = res2.getChain();
    
    if ( chainId1 < chainId2 ) {
      return new Pair<Residue,Residue>(res1,res2);
    }  
    if (chainId2 < chainId1) {
      return new Pair<Residue,Residue>(res2,res1);
    } 
    Integer resNum1 = res1.getResNumber();
    Integer resNum2 = res2.getResNumber();
      
    if ( resNum1 < resNum2 ) {
      return new Pair<Residue,Residue>(res1,res2);
    } else if (resNum2 < resNum1) {
      return new Pair<Residue,Residue>(res2,res1);
    } else {
      return null;
    }
    
  }

}
