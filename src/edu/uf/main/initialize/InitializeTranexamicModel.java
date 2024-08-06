package edu.uf.main.initialize;

import edu.uf.Diffusion.Diffuse;
import edu.uf.compartments.GridFactory;
import edu.uf.compartments.Voxel;
import edu.uf.control.MultiThreadExec;
import edu.uf.interactable.Heme;
import edu.uf.interactable.IL10;
import edu.uf.interactable.Iron;
import edu.uf.interactable.Lactoferrin;
import edu.uf.interactable.MIP1B;
import edu.uf.interactable.MIP2;
import edu.uf.interactable.Macrophage;
import edu.uf.interactable.Neutrophil;
import edu.uf.interactable.PneumocyteII;
import edu.uf.interactable.TGFb;
import edu.uf.interactable.TNFa;
import edu.uf.interactable.Transferrin;
import edu.uf.interactable.Afumigatus.TAFC;
import edu.uf.interactable.Afumigatus.TranexamicAcid;
import edu.uf.utils.Constants;

public class InitializeTranexamicModel extends InitializeBaseModel{

    public void initializeMolecules(Diffuse diffuse, boolean verbose) {
    	if(verbose) {
    		System.out.println("Initializing Iron, TAFC, Lactoferrin, Transferrin, TNF-a, IL10, TGF-b, MIP2, MIP1-b, Heme, and Tranexamic-Acid");
    	}
    	
    	int xbin = GridFactory.getXbin();
    	int ybin = GridFactory.getYbin();
    	int zbin = GridFactory.getZbin();
    	
    	Iron iron = Iron.getMolecule(new double[1][xbin][ybin][zbin], null, new int[] {});
    	TAFC tafc = TAFC.getMolecule(new double[2][xbin][ybin][zbin], diffuse, new int[] {});
    	Lactoferrin lactoferrin = Lactoferrin.getMolecule(new double[3][xbin][ybin][zbin], diffuse, new int[] {Neutrophil.ACTIVE});
    	Transferrin transferrin = Transferrin.getMolecule(new double[3][xbin][ybin][zbin], diffuse, new int[] {});
    	TNFa tnfa = TNFa.getMolecule(new double[1][xbin][ybin][zbin], diffuse, new int[] {Macrophage.M1, Macrophage.M2B, PneumocyteII.ACTIVE, PneumocyteII.MIX_ACTIVE});
    	IL10 il10 = IL10.getMolecule(new double[1][xbin][ybin][zbin], diffuse, new int[] {Macrophage.M1, Macrophage.M2A, Macrophage.M2B, Macrophage.M2C});
    	TGFb tgfb = TGFb.getMolecule(new double[1][xbin][ybin][zbin], diffuse, new int[] {Macrophage.M2C});
    	MIP2 mip2 = MIP2.getMolecule(new double[1][xbin][ybin][zbin], diffuse, new int[] {Macrophage.M1, PneumocyteII.ACTIVE});
    	MIP1B mip1b = MIP1B.getMolecule(new double[1][xbin][ybin][zbin], diffuse, new int[] {Macrophage.M1, PneumocyteII.ACTIVE});
    	Heme heme = Heme.getMolecule(new double[1][xbin][ybin][zbin], null, new int[] {});
    	
    	
    	
    	MultiThreadExec.setMolecule(iron);
    	MultiThreadExec.setMolecule(tafc);
    	MultiThreadExec.setMolecule(lactoferrin);
    	MultiThreadExec.setMolecule(transferrin);
    	MultiThreadExec.setMolecule(tnfa);
    	MultiThreadExec.setMolecule(il10);
    	MultiThreadExec.setMolecule(tgfb);
    	MultiThreadExec.setMolecule(mip2);
    	MultiThreadExec.setMolecule(mip1b);
    	MultiThreadExec.setMolecule(heme);
    	
    	
    	
    	
    	Voxel.setMolecule(Iron.NAME, iron, true, true);
    	Voxel.setMolecule(TAFC.NAME, tafc, true, true);
    	Voxel.setMolecule(Lactoferrin.NAME, lactoferrin, false, true);
    	for(int x = 0; x < xbin; x++) 
        	for(int y = 0; y < ybin; y++)
        		for(int z = 0; z < zbin; z++) {
        			transferrin.set(Constants.DEFAULT_APOTF_CONCENTRATION, 0, x, y, x);
        			transferrin.set(Constants.DEFAULT_TFFE_CONCENTRATION, 1, x, y, x);
        			transferrin.set(Constants.DEFAULT_TFFE2_CONCENTRATION, 2, x, y, x);
        		}
    	Voxel.setMolecule(Transferrin.NAME, transferrin, false, true);
    	Voxel.setMolecule(TNFa.NAME, tnfa);
    	Voxel.setMolecule(IL10.NAME, il10);
    	Voxel.setMolecule(TGFb.NAME, tgfb);
    	Voxel.setMolecule(MIP2.NAME, mip2);
    	Voxel.setMolecule(MIP1B.NAME, mip1b);
    	Voxel.setMolecule(Heme.NAME, heme, true, false);
    	
    	
    	
    	
    	//this.setSecretionPhenotypes();
    }
    
    public void initializeTranexamicAcid(int[] interations) {
    	TranexamicAcid tx = TranexamicAcid.getMolecule(interations);
    	MultiThreadExec.setMolecule(tx);
    	Voxel.setMolecule(TranexamicAcid.NAME, tx);
    }
}
