package org.jiserte.biopdb.contactmap;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import pair.Pair;
import org.jiserte.biopdb.contacts.ContactCriteria;
import org.jiserte.biopdb.structures.Chain;
import org.jiserte.biopdb.structures.Residue;
import org.jiserte.biopdb.structures.SpacePoint;

/**
 * Calculates the contact map of a pdb and returns a list.
 * A ContactCriteria must be given in order to decide if two residues are in 
 * contact.
 * @author javier
 *
 */
public class ContactMapListGenerator {
	////////////////////////////////////////////////////////////////////////////
	// Public Interface
	/**
	 * Retrieves a list of atomic coordinates in contact for a given org.jiserte.biopdb.
	 * @param pdb
	 * @param criteria
	 * @return
	 */
	public List<Pair<SpacePoint,SpacePoint>> getSpacePointsInContact(
	    Map<Character,Chain> pdb, ContactCriteria criteria) {
		
		List<Pair<SpacePoint,SpacePoint>> result = 
		    new ArrayList<Pair<SpacePoint,SpacePoint>>();
		
		if (criteria.canProvideCandidates()) {
		  
		  Set<Pair<Residue, Residue>> candidates = criteria.getCandidates();
		  
		  for (Pair<Residue, Residue> candidate : candidates) {
		    
        Pair<SpacePoint, SpacePoint> pair = criteria.areInContact(
            candidate.getFirst(), candidate.getSecond());
        
        if (pair!=null) {
          
          result.add(pair);
          
        }
		  }
		  
		} else {

  		List<Residue> residues = new ArrayList<Residue>();
  		
  		for (Chain chain : pdb.values()) {
  
  			residues.addAll(chain.getResiduesCollection());
  			
  		}
  			
  		for (int i =0 ; i< residues.size() -1 ; i++) {
  			
  			for (int j =i+1 ;j< residues.size() ; j++) {
  
  				Pair<SpacePoint, SpacePoint> pair = criteria.areInContact(
  				    residues.get(i), residues.get(j));
  				
  				if (pair!=null) {
  					
  					result.add(pair);
  					
  				}
  				
  			}
  			
  		}
		}
		return result;
		
	}
	////////////////////////////////////////////////////////////////////////////

}
