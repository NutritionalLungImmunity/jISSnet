package edu.uf.utils;

public class Id {
	private static int ID = 0;
	private static int moleculeID = 0;
	
	public static int getId() {
		Id.ID = Id.ID + 1;
		return Id.ID; 
	}
	
	public static int getMoleculeId() {
		Id.moleculeID = Id.moleculeID + 1;
		return Id.moleculeID; 
	}
}
