package edu.uf.intracellularState;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import edu.uf.interactable.MIP2;
import edu.uf.interactable.Neutrophil;
import edu.uf.interactable.TLRBinder;
import edu.uf.interactable.Afumigatus.Afumigatus;
import edu.uf.utils.Constants;
import edu.uf.utils.Rand;

public class NeutrophilStateModel extends BooleanNetwork{

	public static final int size = 4;
	//public static final int NUM_RECEPTORS = 4;
	
	public static final int ACTIVE = 0;
	public static final int DYING = 1;
	public static final int APOPTOTIC = 2;
	public static final int NETOTIC = 3;
	
	{
		this.inputs = new int[NUM_RECEPTORS];
		this.booleanNetwork = new int[size];
	}
	
	
	@Override
	public void processBooleanNetwork(int... args) {
		
		if(this.booleanNetwork[APOPTOTIC] > 0 || this.booleanNetwork[NETOTIC] == 1) {
			//System.out.println(this.booleanNetwork[APOPTOTIC] + " " + this.booleanNetwork[NETOTIC] + " " + this.getPhenotype().get((NeutrophilStateModel.NETOTIC)));
			return;
		}

		int k = 0;
		List<Integer> array = new ArrayList<>(size);
		for(int i = 0; i < size; i++)
			array.add(i);
		while(true) {
			if(k++ > 3)break;
			Collections.shuffle(array, new Random());
			for(int i : array) {
				switch(i) {
					case 0:
						this.booleanNetwork[ACTIVE] = input(Afumigatus.DEF_OBJ) | input(TLRBinder.getBinder()) | this.booleanNetwork[ACTIVE];
						break;
					case 1:
						
						this.booleanNetwork[DYING] = this.booleanNetwork[ACTIVE] > 0 ? (Rand.getRand().randunif() < Constants.NEUTROPHIL_HALF_LIFE ? 1 : 0) : 0;
						break;
					case 2:
						if(this.booleanNetwork[APOPTOTIC] == 0) {
							int j = this.booleanNetwork[DYING] > 0 ? (Rand.getRand().randunif() < Constants.PR_NEUT_APOPT ? 2 : 1) : 0;
						
							this.booleanNetwork[APOPTOTIC] =  j == 2 ? 1 : 0;
							this.booleanNetwork[NETOTIC] = j == 1 ? 1 : 0; //
							//if(this.booleanNetwork[DYING] > 0) System.out.println(this.booleanNetwork[NETOTIC]  + " " + this.booleanNetwork[APOPTOTIC]);
						}
						break;
					case 3:
						this.booleanNetwork[APOPTOTIC] = this.booleanNetwork[ACTIVE] == 0 ? (Rand.getRand().randunif() < Constants.NEUTROPHIL_HALF_LIFE ? 1 : 0) : 0;
						//System.out.println(this.booleanNetwork[APOPTOTIC]);
						break;
					default:
						System.err.println("No such interaction " + i + "!");
						break;
				}
			}
		}
		
		for(int i = 0; i < NUM_RECEPTORS; i++)
			this.inputs[i] = 0;
		
		this.clearPhenotype();
		
		if(this.booleanNetwork[ACTIVE] > 0) {
			this.getPhenotype().put(NeutrophilStateModel.ACTIVE, this.booleanNetwork[ACTIVE]);
		} 
		if(this.booleanNetwork[APOPTOTIC] == 1) {
			this.getPhenotype().put(NeutrophilStateModel.APOPTOTIC, this.booleanNetwork[APOPTOTIC]);
		}
		if(this.booleanNetwork[NETOTIC] == 1) {
			this.getPhenotype().put(NeutrophilStateModel.NETOTIC, this.booleanNetwork[NETOTIC]);
		}
			
	
		
	}
	
}
