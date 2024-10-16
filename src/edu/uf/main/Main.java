package edu.uf.main;

import java.io.File;

import edu.uf.Diffusion.Diffuse;
import edu.uf.Diffusion.FADIPeriodic;
import edu.uf.compartments.GridFactory;
import edu.uf.compartments.MacrophageRecruiter;
import edu.uf.compartments.NeutrophilRecruiter;
import edu.uf.compartments.Recruiter;
import edu.uf.interactable.Neutrophil;
import edu.uf.interactable.Afumigatus.Afumigatus;
import edu.uf.main.initialize.Initialize;
import edu.uf.main.initialize.InitializeBaseModel;
import edu.uf.main.print.PrintHemeModel;
import edu.uf.main.run.Run;
import edu.uf.main.run.RunSingleThread;
import edu.uf.utils.Constants;

public class Main {
	
	private static void baseModel(String[] args) throws InterruptedException {
		

		Neutrophil.MECH = Neutrophil.KILL_INJ; //use this for kill-injured (NET only kills type-I pneumocytes previously injured by hyphae)
		//Neutrophil.MECH = Neutrophil.KILL; //use this for kill (NET kill type-I pneumocytes with a probability independent of hyphae)
		//Neutrophil.MECH = Neutrophil.DNASE; //use this for no NETs (NET will not kill type-I pneumocytes)
		//the three above are mutually exclusive. "KILL_INJ" is the default.

		
		Constants.PR_NET_KILL_EPI *= 0.25;//Double.parseDouble(args[0]);//(0.1*Double.parseDouble(args[0]));
		Constants.NET_COUNTER_INHIBITION = 0.0;//Double.parseDouble(args[0]);
		Constants.HEME_QTTY *= 10;//Double.parseDouble(args[0]);
		Constants.HEME_UP  *= 0.75;//0.75;//Double.parseDouble(args[1]);
		
		int i = 0;//Integer.parseInt(args[0]);
		
		String filename = "HemeNET025_";// + args[0] + "_";

		filename += i + ".tsv";
		 	
		Initialize initialize = new InitializeBaseModel();
		Run run = new RunSingleThread();
		PrintHemeModel stat = new PrintHemeModel();
		
		
		int xbin = 10;
		int ybin = 10;
		int zbin = 10;
        
        //int pne = (int) (xbin * ybin * zbin  * 0.64);
        
        String[] input = new String[]{"0", "1920", "15", "640"};
        
        //Constants.Kd_Granule = 1e12;
        
        /*Constants.ITER_TO_GROW = 328/2;//Integer.parseInt(args[0]);  //REF 328
        Constants.PR_N_HYPHAE = 0;//Double.parseDouble(args[1]) * Constants.PR_N_HYPHAE ;
        Constants.PR_MA_HYPHAE = Constants.PR_MA_HYPHAE/4;//Double.parseDouble(args[1]) * Constants.PR_MA_HYPHAE;
        Constants.GRANULE_QTTY = 1;//Double.parseDouble(args[2]);
        Constants.Kd_Granule = 0.1 * 1e12;//10 * 1e12; 3.98 //REF:10
        //Constants.Kd_Granule = Double.parseDouble(args[0]) * 1e12;
        Constants.Granule_HALF_LIFE = -1+Math.log(0.5)/20;//1+Math.log(0.5)/20.0;// //harmonic avg of the best //REF:20
        Constants.MAX_MA = 660;//Integer.parseInt(args[3]);
        Constants.MAX_N = 660;//Integer.parseInt(args[4]);
        
        Constants.PR_DEPLETION = 0.0;*/
        
        //Constants.P_TNF_QTTY = 1*Constants.P_TNF_QTTY;
        //Constants.MIN_MA = 15;
        
        //heme diffusion is off!!!
        
        //Constants.PR_ASPERGILLUS_CHANGE = 0.1732868;
        Constants.MAX_N = 522;//*10 for WT
        //Constants.MAX_MA = 0;//0.75;//*0.75;
        Constants.RECRUITMENT_RATE_N *= 1.0;//0.25;//0.125;//0.25;
        Constants.NEUTROPHIL_HALF_LIFE = 0.05776227;
        Constants.NET_HALF_LIFE = 0.007701635;
        Constants.DNAse_KCAT = 0.007701635 * 10.0;
        
        Constants.DNAse_HALF_LIFE = 0.999;
        
        //Constants.ITER_TO_GROW =  30;
        
        Constants.PR_ASP_KILL_EPI *= 1.0;
        
        double f = 0.1;
        double pdeFactor = Constants.D/(30/Constants.TIME_STEP_SIZE); 
        double dt = 1;
        Diffuse diffusion = new FADIPeriodic(f, pdeFactor, dt);
        
        initialize.createPeriodicGrid(xbin, ybin, zbin);
        initialize.initializeMolecules(diffusion, false);
        initialize.initializePneumocytes(Integer.parseInt(input[3]));
        //initialize.initializeLiver(grid, xbin, ybin, zbin);
        initialize.initializeMacrophage(Integer.parseInt(input[2]));
        initialize.initializeNeutrophils(0);
        initialize.initializeTypeIPneumocytes(Integer.parseInt(input[3])/2);
        initialize.initializeBlood();
        initialize.infect(Integer.parseInt(input[1]), Afumigatus.RESTING_CONIDIA, Constants.CONIDIA_INIT_IRON, -1, false);
        stat.grid = GridFactory.getGrid();

        Recruiter[] recruiters = new Recruiter[2];
        recruiters[0] = new MacrophageRecruiter();
        recruiters[1] = new NeutrophilRecruiter();
        //recruiters[0] = new MacrophageReplenisher();
        //recruiters[1] = new NeutrophilReplenisher();
        
        run.run(
        		1441,//(int) 72*2*(30/Constants.TIME_STEP_SIZE),  
        		xbin, 
        		ybin, 
        		zbin,  
        		recruiters,
        		false,
        		null, //new File(filename),
        		-1,
        		stat
        );
        
        stat.close();
        //System.out.println((toc - tic));
	}

	public static void main(String[] args) throws Exception {
		//Turnover and degrade changed: Turnover now is automatic. In the future there should be turnover and then degrade in the serum.
		//Current Turnover rate is wrong. Should it be (1-Turnover_rate)?
		System.out.println("jISSnet");
		long tic = System.currentTimeMillis();
		Main.baseModel(args);
		long toc = System.currentTimeMillis();
		System.out.println((toc - tic));
	}

}


//15310, 28325, 54318