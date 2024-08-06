/*package edu.uf.interactable.invitro;

import edu.uf.Diffusion.Diffuse;
import edu.uf.interactable.Afumigatus;
import edu.uf.interactable.Interactable;
import edu.uf.interactable.Molecule;
import edu.uf.intracellularState.Molecules;
import edu.uf.utils.Constants;
import edu.uf.utils.Rand;
import edu.uf.utils.Util;

public class Glutamate extends Molecule implements Invitro{

	public static final String NAME = Molecules.Glutamate;
	public static final int NUM_STATES = 1;
	
	private static Glutamate molecule = null; 
    
    private Glutamate(double[][][][] qttys, Diffuse diffuse) {
		super(qttys, diffuse);
	}
    
    public static Glutamate getMolecule(double[][][][] values, Diffuse diffuse) {
    	if(molecule == null) {
    		molecule = new Glutamate(values, diffuse);
    	}
    	return molecule;
    }
    
    public static Glutamate getMolecule() {
    	return molecule;
    }
    

    public int getIndex(String str) {
        return 0;
    }

    public void degrade(int x, int y, int z) {}

    public void computeTotalMolecule(int x, int y, int z) {
    	this.totalMoleculesAux[0] = this.totalMoleculesAux[0] + this.values[0][x][y][z];
    }

    protected boolean templateInteract(Interactable interactable, int x, int y, int z) {
        if (interactable instanceof Afumigatus) {
            InvitroAfumigatus afumigatus = (InvitroAfumigatus) interactable;
            if(afumigatus.getStatus() != Afumigatus.RESTING_CONIDIA) {
            	if(Util.activationFunction(this.get(x, y, z), Constants.Kd_GLU, Constants.STD_UNIT_T) > Rand.getRand().randunif()) {
            		double qtty = Constants.GLU_UPTAKE_RATE * this.get(x, y, z);
            		this.dec(qtty, x, y, z);
            		afumigatus.setBooleanNetwork(InvitroAfumigatus.Glu, 1);
            	}
            }
            return true;  
        }
        return interactable.interact(this, x, y, z);
    }

	@Override
	public String getName() {
		return NAME;
	}

	@Override
	public double getThreshold() {
		return 0;
	}

	@Override
	public int getNumState() {
		return NUM_STATES;
	}
	
	@Override
	public String getReceptor() {
		return null; //REVIEW
	}
	
}*/
