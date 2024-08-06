package edu.uf.main.print;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Map.Entry;

import edu.uf.compartments.Voxel;
import edu.uf.interactable.IL10;
import edu.uf.interactable.IL6;
import edu.uf.interactable.Interactable;
import edu.uf.interactable.MIP2;
import edu.uf.interactable.Macrophage;
import edu.uf.interactable.Neutrophil;
import edu.uf.interactable.PneumocyteII;
import edu.uf.interactable.TGFb;
import edu.uf.interactable.TNFa;
import edu.uf.interactable.Afumigatus.Afumigatus;

public class PrintBaseModel extends PrintStat{
	
	int PResting;
	int PMixActive;
	int PActive;
	int MAResting;
	int MAMixActive;
	int MAActive;
	int MAInactive;
	
	public Voxel[][][] grid;

	@Override
	public void printStatistics(int k, File file){
		count();
		if(k%15 != 0)return;
		String str = k + "\t" + 
	              Afumigatus.getTotalCells() + "\t" +
	              Afumigatus.getTotalRestingConidia() + "\t" +
	              Afumigatus.getTotalSwellingConidia() + "\t" +
	              Afumigatus.getTotalGerminatingConidia() + "\t" +
	              Afumigatus.getTotalHyphae() + "\t" +
	              //PResting + "\t" +
	              PMixActive + "\t" +
	              PActive + "\t" +
	              /*(TAFC.getMolecule().getTotalMolecule(0) + TAFC.getMolecule().getTotalMolecule(1)) + "\t" +
	              TAFC.getMolecule().getTotalMolecule(0) + "\t" +
	              TAFC.getMolecule().getTotalMolecule(1) + "\t" +
	              Lactoferrin.getMolecule().getTotalMolecule(0) + "\t" +
	              Lactoferrin.getMolecule().getTotalMolecule(1) + "\t" +
	              Lactoferrin.getMolecule().getTotalMolecule(2) + "\t" +
	              (Transferrin.getMolecule().getTotalMolecule(0) + Transferrin.getMolecule().getTotalMolecule(1) + Transferrin.getMolecule().getTotalMolecule(2)) + "\t" +
	              Transferrin.getMolecule().getTotalMolecule(0) + "\t" +
	              Transferrin.getMolecule().getTotalMolecule(1) + "\t" +
	              Transferrin.getMolecule().getTotalMolecule(2) + "\t" +
	              Hepcidin.getMolecule().getTotalMolecule(0) + "\t" +*/
	              TGFb.getMolecule().getTotalMolecule(0) + "\t" +
	              IL6.getMolecule().getTotalMolecule(0) + "\t" +
	              IL10.getMolecule().getTotalMolecule(0) + "\t" +
	              TNFa.getMolecule().getTotalMolecule(0) + "\t" +
	              //MIP1B.getMolecule().getTotalMolecule(0) + "\t" +
	              MIP2.getMolecule().getTotalMolecule(0) + "\t" +
	              //Erythrocyte.getTotalCells() + "\t" + 
	              /*MAResting + "\t" +
	              MAMixActive + "\t" +
	              MAActive + "\t" +
	              MAInactive + "\t" +*/
	              Macrophage.getTotalCells() + "\t" +
	              Neutrophil.getTotalCells();
		
		if(file == null) {
			System.out.println(str);
		}else {
			try {
		
				if(getPrintWriter() == null) 
					setPrintWriter(new PrintWriter(file)); 
				getPrintWriter().println(str);
			}catch(FileNotFoundException ex) {
				ex.printStackTrace();
			}
		}
	}
	
	void count() {
		PResting = 0;
		PMixActive = 0;
		PActive = 0;
		MAResting = 0;
		MAMixActive = 0;
		MAActive = 0;
		MAInactive = 0;
		for(Voxel[][] VV : grid) {
			for(Voxel[] V : VV) {
				for(Voxel v : V) {
					for(Entry<Integer, Interactable> entry : v.getInteractables().entrySet()) {
						Interactable cell = entry.getValue();
						if(cell instanceof PneumocyteII) {
							PneumocyteII p = (PneumocyteII) cell;
							if(p.hasPhenotype(PneumocyteII.MIX_ACTIVE))PMixActive++;
							else if(p.hasPhenotype(PneumocyteII.ACTIVE))PActive++;
						}/*else if(cell instanceof Macrophage) {
							Macrophage p = (Macrophage) cell;
							if(p.getPhenotype() == Phenotypes.RESTING)MAResting++;
							else if(p.getPhenotype() == Phenotypes.MIX_ACTIVE)MAMixActive++;
							else if(p.getPhenotype() == Phenotypes.ACTIVE)MAActive++;
							else if(p.getPhenotype() == Phenotypes.INACTIVE)MAInactive++;
						}*/
					}
				}
			}
		}
	}
	

}


//307200
