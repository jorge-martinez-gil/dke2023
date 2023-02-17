
package dke2023;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.moeaframework.Executor;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.core.variable.RealVariable;
import org.moeaframework.problem.AbstractProblem;

import com.github.chen0040.data.utils.TupleTwo;
import com.github.chen0040.gp.commons.BasicObservation;
import com.github.chen0040.gp.commons.Observation;

import jmetal.core.Variable;
import net.sourceforge.jFuzzyLogic.FIS;

public class dke2023 {

public static class dke extends AbstractProblem {
		

		
		static double [] source_test;
		static double [] target_test;
		static double [] source_training;
		static double [] target_training;
		
	    private static String [] valueS = {"null", "poor", "good", "excellent"};
	    private static String [] valueAS = {"COG", "MM", "LM", "RM"};
		
		static String [] v = null;
		int i = 0;
		
		static List<Observation> trainingData;
		public static List<Observation> testingData;

		/**
		 * Constructs a new instance of the program
		 */
		public dke() {
			super(62, 3);
			
			List<Observation> data = testingData = generate("datasets/mc.txt");
			CollectionUtils.shuffle(data);
			TupleTwo<List<Observation>, List<Observation>> split_data = CollectionUtils.split(data, 0.8);
			trainingData = split_data._1();
			testingData = split_data._2();
			
		}
		
		/**
		 * Load the data from the neural network in main memory
		 */
		private static List<Observation> generate(String file){
			      
			List<Observation> result = new ArrayList<>();

			        try {
			            final FileInputStream stream = new FileInputStream(file);
						try (BufferedReader reader = new BufferedReader(new InputStreamReader(stream))) {
							String data = "";
							do {
							    data = reader.readLine();
							    if (data == null || data.trim().length() == 0) {
							    	continue;
							    } else {
							    	final String [] line = data.split(",");
							        if (line.length == 8) {
							            try {
							                final double a = Double.parseDouble(line[0]);
							                final double b = Double.parseDouble(line[1]);
							                final double c = Double.parseDouble(line[2]);
							                final double d = Double.parseDouble(line[3]);
							                final double e = Double.parseDouble(line[4]);
							                

								            Observation observation = new BasicObservation(4, 1);
								            observation.setInput(0, b);
								            observation.setInput(1, c);
								            observation.setInput(2, d);
								            observation.setInput(3, e);
								            observation.setOutput(0, a);
								            result.add(observation);
							                
							                    
							            } catch(final NumberFormatException e) {
							                System.out.println(e);
							            }
							        }
							    } 
							} while(data != null);
						}            

			            
			        } catch(Exception e) {
			            e.printStackTrace();
			        }
			        
			        
			      return result;
			   }
		
		
		/** 
		 * Compute the cosineSimilarity of two vectors
		 */ 
		public static double cosineSimilarity(double[] vectorA, double[] vectorB) {
		    double dotProduct = 0.0;
		    double normA = 0.0;
		    double normB = 0.0;
		    for (int i = 0; i < vectorA.length; i++) {
		        dotProduct += vectorA[i] * vectorB[i];
		        normA += Math.pow(vectorA[i], 2);
		        normB += Math.pow(vectorB[i], 2);
		    }   
		    return dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
		}
		
		
		/** 
		 * Compute the Pearson Correlation Coefficient
		 */ 
		public static double getPearson(double[] scores1, double[] scores2) {
			double result = 0;
			double sum_sq_x = 0;
			double sum_sq_y = 0;
			double sum_coproduct = 0;
			double mean_x = scores1[0];
			double mean_y = scores2[0];
			for (int i = 2; i < scores1.length + 1; i += 1) {
				double sweep = Double.valueOf(i - 1) / i;
				double delta_x = scores1[i - 1] - mean_x;
				double delta_y = scores2[i - 1] - mean_y;
				sum_sq_x += delta_x * delta_x * sweep;
				sum_sq_y += delta_y * delta_y * sweep;
				sum_coproduct += delta_x * delta_y * sweep;
				mean_x += delta_x / i;
				mean_y += delta_y / i;
			}
			double pop_sd_x = (double) Math.sqrt(sum_sq_x / scores1.length);
			double pop_sd_y = (double) Math.sqrt(sum_sq_y / scores1.length);
			double cov_x_y = sum_coproduct / scores1.length;
			result = cov_x_y / (pop_sd_x * pop_sd_y);
			
			if (Double.isNaN(result))
				return -99;
			else
				return result;
		}
		
		@SuppressWarnings("null")
		/** 
		 * Compute the Spearman Rank Correlation Coefficient
		 */ 
		private static double getSpearman(double [] X, double [] Y) {
		       
			try {
			/* Error check */
		        if (X == null || Y == null || X.length != Y.length) {
		            return (Double) null;
		        }
		        
		        /* Create Rank arrays */
		        int [] rankX = getRanks(X);
		        int [] rankY = getRanks(Y);

		        /* Apply Spearman's formula */
		        int n = X.length;
		        double numerator = 0;
		        for (int i = 0; i < n; i++) {
		            numerator += Math.pow((rankX[i] - rankY[i]), 2);
		        }
		        numerator *= 6;
		        return 1 - numerator / (n * ((n * n) - 1));
			}
			catch (Exception e) {
				return -999;
			}
		    }
		    
		    /* Returns a new (parallel) array of ranks. Assumes unique array values */
		    public static int[] getRanks(double [] array) {
		        int n = array.length;
		        
		        /* Create Pair[] and sort by values */
		        Pair [] pair = new Pair[n];
		        for (int i = 0; i < n; i++) {
		            pair[i] = new Pair(i, array[i]);
		        }
		        Arrays.sort(pair, new PairValueComparator());

		        /* Create and return ranks[] */
		        int [] ranks = new int[n];
		        int rank = 1;
		        for (Pair p : pair) {
		            ranks[p.index] = rank++;
		        }
		        return ranks;
		    }
		
		/** 
		* Construct the fuzzy program
		*/ 
		public static void replace(String [] target){  
			    try {  
			        FileReader fr = new FileReader("fcl/tipperB.fcl");  
			        BufferedReader br = new BufferedReader(fr);  
			        FileWriter fw = new FileWriter("fcl/tipper99.fcl");  
			        PrintWriter bw = new PrintWriter(fw);  
			        String line = null;  
			  
			        while((line=br.readLine()) != null) {  
			        	line = line.replaceAll("aaa", String.valueOf(target[0]));    
			        	line = line.replaceAll("bbb", String.valueOf(target[1])); 
			        	line = line.replaceAll("ccc", String.valueOf(target[2])); 
			        	line = line.replaceAll("ddd", String.valueOf(target[3])); 
			        	line = line.replaceAll("eee", String.valueOf(target[4])); 
			        	line = line.replaceAll("fff", String.valueOf(target[5])); 
			        	line = line.replaceAll("ggg", String.valueOf(target[6])); 
			        	line = line.replaceAll("hhh", String.valueOf(target[7])); 
			        	line = line.replaceAll("iii", String.valueOf(target[8])); 
			        	line = line.replaceAll("jjj", String.valueOf(target[9])); 
			        	line = line.replaceAll("kkk", String.valueOf(target[10]));
			        	line = line.replaceAll("lll", String.valueOf(target[11])); 
			        	line = line.replaceAll("mmm", String.valueOf(target[12])); 
			        	line = line.replaceAll("nnn", String.valueOf(target[13])); 
			        	line = line.replaceAll("ooo", String.valueOf(target[14])); 
			        	line = line.replaceAll("ppp", String.valueOf(target[15])); 
			        	line = line.replaceAll("qqq", String.valueOf(target[16])); 
			        	line = line.replaceAll("rrr", String.valueOf(target[17])); 
			        	
			        	line = line.replaceAll("2aa2", String.valueOf(target[18]));    
			        	line = line.replaceAll("2bb2", String.valueOf(target[19])); 
			        	line = line.replaceAll("2cc2", String.valueOf(target[20])); 
			        	line = line.replaceAll("2dd2", String.valueOf(target[21])); 
			        	line = line.replaceAll("2ee2", String.valueOf(target[22])); 
			        	line = line.replaceAll("2ff2", String.valueOf(target[23])); 
			        	line = line.replaceAll("2gg2", String.valueOf(target[24])); 
			        	line = line.replaceAll("2hh2", String.valueOf(target[25])); 
			        	line = line.replaceAll("2ii2", String.valueOf(target[26])); 
			        	line = line.replaceAll("2jj2", String.valueOf(target[27])); 
			        	line = line.replaceAll("2kk2", String.valueOf(target[28]));
			        	line = line.replaceAll("2ll2", String.valueOf(target[29])); 
			        	line = line.replaceAll("2mm2", String.valueOf(target[30])); 
			        	line = line.replaceAll("2nn2", String.valueOf(target[31])); 
			        	line = line.replaceAll("2oo2", String.valueOf(target[32])); 
			        	line = line.replaceAll("2pp2", String.valueOf(target[33])); 
			        	line = line.replaceAll("2qq2", String.valueOf(target[34])); 
			        	line = line.replaceAll("2rr2", String.valueOf(target[35])); 
			        	
			        	line = line.replaceAll("3aa3", String.valueOf(target[36]));    
			        	line = line.replaceAll("3bb3", String.valueOf(target[37])); 
			        	line = line.replaceAll("3cc3", String.valueOf(target[38])); 
			        	line = line.replaceAll("3dd3", String.valueOf(target[39])); 
			        	line = line.replaceAll("3ee3", String.valueOf(target[40])); 
			        	line = line.replaceAll("3ff3", String.valueOf(target[41])); 
			        	line = line.replaceAll("3gg3", String.valueOf(target[42])); 
			        	line = line.replaceAll("3hh3", String.valueOf(target[43])); 
			        	line = line.replaceAll("3ii3", String.valueOf(target[44])); 
			        	line = line.replaceAll("3jj3", String.valueOf(target[45])); 
			        	line = line.replaceAll("3kk3", String.valueOf(target[46]));
			        	line = line.replaceAll("3ll3", String.valueOf(target[47])); 
			        	line = line.replaceAll("3mm3", String.valueOf(target[48])); 
			        	line = line.replaceAll("3nn3", String.valueOf(target[49])); 
			        	line = line.replaceAll("3oo3", String.valueOf(target[50])); 
			        	line = line.replaceAll("3pp3", String.valueOf(target[51])); 
			        	line = line.replaceAll("3qq3", String.valueOf(target[52])); 
			        	line = line.replaceAll("3rr3", String.valueOf(target[53])); 
			        	
			        	line = line.replaceAll("4aa4", String.valueOf(target[54]));    
			        	line = line.replaceAll("4bb4", String.valueOf(target[55])); 
			        	line = line.replaceAll("4cc4", String.valueOf(target[56])); 
			        	line = line.replaceAll("4dd4", String.valueOf(target[57])); 
			        	line = line.replaceAll("4ee4", String.valueOf(target[58])); 
			        	line = line.replaceAll("4ff4", String.valueOf(target[59])); 
			        	line = line.replaceAll("4gg4", String.valueOf(target[60])); 
			        	line = line.replaceAll("4hh4", String.valueOf(target[61])); 
			        	
			        
			        	bw.println(line);  
			            } 
			        
			        fr.close();
			        fw.close();
			        
			      
			        }  
			    catch(Exception e){e.printStackTrace();} 
		   }
		   
		/**
		 * Copies the FLC
		 */
		private static void copyFileUsingStream(String s, String d) throws IOException {
			    
	        File inputFile = new File(s);
	        File tempFile = new File("fcl/temp");

	        BufferedReader reader = new BufferedReader(new FileReader(inputFile));
	        BufferedWriter writer = new BufferedWriter(new FileWriter(tempFile));

	        String currentLine;

	        while((currentLine = reader.readLine()) != null) {
	            if(currentLine.contains("null") == false)
	            	writer.write(currentLine + System.getProperty("line.separator"));
	        }
	        
	        writer.close(); 
	        reader.close(); 

			File source = new File ("fcl/temp");
			File dest = new File (d);
				
			InputStream is = null;
		    OutputStream os = null;
		    try {
		        is = new FileInputStream(source);
		        os = new FileOutputStream(dest);
		        byte[] buffer = new byte[1024];
		        int length;
		        while ((length = is.read(buffer)) > 0) {
		            os.write(buffer, 0, length);
		        }
		    } finally {
		        is.close();
		        os.close();
		    }
		}
		
		/**
		 * Calculates Validation
		 */
		public static void finalValidation (List<Observation> lista, Solution solution) {
			
			double[] source;
			double[] target;
			   	
			double[] x = EncodingUtils.getReal(solution);
			
			// Complete individual
	        String[] geneString = new String[62];
	        
	        // Rules
	        for (int i = 0; i < 54; i++) {
	        	int num = (int) x[i];
	        	geneString[i] = valueS[num];
	        	if (i == 12)
	        		geneString[i] = valueAS[num];
	        }
	        
	        
	        // Coordinates
	        for (int i = 54; i < 62; i++) {      	
	        	double val = (double)Math.round(x[i]/2.99 * 1000d) / 1000d;
	        	val = val * 0.15;
	        	
	 		    double a = Double.valueOf(val);
	 		    double b = Double.valueOf(val + a);
	 		    
	 		    double c = Double.valueOf(val * 2);
	 		    double d = Double.valueOf(val + c);
	 		    double e = Double.valueOf(val + d);
	 		    double f1 = Double.valueOf( val + e);
	 		    
	 		    double g = Double.valueOf(val * 4);
	 		    double h = Double.valueOf(val + g);
	 		    
	 		    
	 		    geneString[54] = String.valueOf((double)Math.round(a * 1000d) / 1000d);
	 		    geneString[55] = String.valueOf((double)Math.round(b * 1000d) / 1000d);
			    
	 		    geneString[56] = String.valueOf((double)Math.round(c * 1000d) / 1000d);
	 		    geneString[57] = String.valueOf((double)Math.round(d * 1000d) / 1000d);
	 		    geneString[58] = String.valueOf((double)Math.round(e * 1000d) / 1000d);
	 		    geneString[59] = String.valueOf((double)Math.round(f1 * 1000d) / 1000d);
			    
	 		    geneString[60] = String.valueOf((double)Math.round(g * 1000d) / 1000d);
	 		    geneString[61] = String.valueOf((double)Math.round(h * 1000d) / 1000d);
	        
	        }
	        
	        String fileName = "fcl/best.fcl";
			replace (geneString);								 									   	
			FIS fis = new FIS ();
	        try {
	        	copyFileUsingStream("fcl/tipper99.fcl", fileName);
	        	fis = FIS.load(fileName, true);
			} catch (Exception e) {
				return;
			}
			
			int i = 0;
			source = new double [lista.size()];
			target = new double [lista.size()];
			for(Observation observation : lista) {
	 
		        fis.setVariable("alga", observation.getInput(0));
		        fis.setVariable("algb", observation.getInput(1));
		        fis.setVariable("algc", observation.getInput(2));
		        fis.setVariable("algd", observation.getInput(3));
		        fis.evaluate();
		
				double predicted = fis.getVariable("score").defuzzify();
				double actual = observation.getOutput(0);
			 
				source[i] = predicted;
				target[i] = actual;
				i++;
				
				System.out.println (predicted + " " + actual);
			}
			
			System.out.println ("Correl: " + getPearson(target,source));
			
		}
		

		/**
		 * Constructs a new solution and defines the bounds of the decision
		 * variables.
		 */
		@Override
		public Solution newSolution() {
			Solution solution = new Solution(getNumberOfVariables(), 
					getNumberOfObjectives());

			for (int i = 0; i < getNumberOfVariables(); i++) {
				solution.setVariable(i, new RealVariable(0, 3.99));
			}

			return solution;
		}
		

		/**
		 * Evaluation function
		 */
		@Override
		public void evaluate(Solution solution) {
			
			// From 0 to 54 --> Rules
			// From 54 to 62 --> Coordinates
			
			double[] source;
			double[] target;
			
			double[] x = EncodingUtils.getReal(solution);
			double[] f = new double[numberOfObjectives];
			
			i++;
			
			// Complete individual
	        String[] geneString = new String[62];
	        
	        // Rules
	        for (int i = 0; i < 54; i++) {
	        	int num = (int) x[i];
	        	geneString[i] = valueS[num];
	        	if (i == 12)
	        		geneString[i] = valueAS[num];
	        }
	        
	        
	        // Coordinates
	        for (int i = 54; i < 62; i++) {      	
	        	double val = (double)Math.round(x[i]/2.99 * 1000d) / 1000d;
	        	val = val * 0.15;
	        	
	 		    double a = Double.valueOf(val);
	 		    double b = Double.valueOf(val + a);
	 		    
	 		    double c = Double.valueOf(val * 2);
	 		    double d = Double.valueOf(val + c);
	 		    double e = Double.valueOf(val + d);
	 		    double f1 = Double.valueOf( val + e);
	 		    
	 		    double g = Double.valueOf(val * 4);
	 		    double h = Double.valueOf(val + g);
	 		    
	 		    
	 		    geneString[54] = String.valueOf((double)Math.round(a * 1000d) / 1000d);
	 		    geneString[55] = String.valueOf((double)Math.round(b * 1000d) / 1000d);
			    
	 		    geneString[56] = String.valueOf((double)Math.round(c * 1000d) / 1000d);
	 		    geneString[57] = String.valueOf((double)Math.round(d * 1000d) / 1000d);
	 		    geneString[58] = String.valueOf((double)Math.round(e * 1000d) / 1000d);
	 		    geneString[59] = String.valueOf((double)Math.round(f1 * 1000d) / 1000d);
			    
	 		    geneString[60] = String.valueOf((double)Math.round(g * 1000d) / 1000d);
	 		    geneString[61] = String.valueOf((double)Math.round(h * 1000d) / 1000d);
	        
	        }
	        
	        String fileName = "fcl/result.fcl";
			replace (geneString);								 									   	
			FIS fis = new FIS ();
	        try {
	        	copyFileUsingStream("fcl/tipper99.fcl", fileName);
	        	fis = FIS.load(fileName, true);
			} catch (Exception e) {
				return;
			}

	        
			int i = 0;
			source = new double [trainingData.size()];
			target = new double [trainingData.size()];
			for(Observation observation : trainingData) {
	 
		        fis.setVariable("alga", observation.getInput(0));
		        fis.setVariable("algb", observation.getInput(1));
		        fis.setVariable("algc", observation.getInput(2));
		        fis.setVariable("algd", observation.getInput(3));
		        fis.evaluate();
		
				double predicted = fis.getVariable("score").defuzzify();
				double actual = observation.getOutput(0);

				//System.out.println (predicted + " " + actual);
				 
				source[i] = predicted;
				target[i] = actual;
				i++;
			}

	        f[0] = -getPearson(source, target);
	        //f[2] = cosineSimilarity (source, target);
	        
			i = 0;
			source = new double [testingData.size()];
			target = new double [testingData.size()];
			for(Observation observation : testingData) {
	 
		        fis.setVariable("alga", observation.getInput(0));
		        fis.setVariable("algb", observation.getInput(1));
		        fis.setVariable("algc", observation.getInput(2));
		        fis.setVariable("algd", observation.getInput(3));
		        fis.evaluate();
		
				double predicted = fis.getVariable("score").defuzzify();
				double actual = observation.getOutput(0);

				//System.out.println (predicted + " " + actual);
				 
				source[i] = predicted;
				target[i] = actual;
				i++;
			}

	        //f[1] = -getPearson(source, target);
	        
	        solution.setObjectives(f);
		}
	}


	//main
	public static void main(String[] args) {
		//configure and run the main function
		long startTime = System.nanoTime();
		NondominatedPopulation result = new Executor()
				.withProblemClass(dke.class)
				.withAlgorithm("DE")
				.withMaxEvaluations(500)
				.run();
		 
		long endTime = System.nanoTime();
		long timeElapsed = endTime - startTime;
		System.out.println("Execution time in milliseconds : " + timeElapsed / 1000000);

		System.out.format("Training		Test \n");
		
		for (Solution solution : result) {
			System.out.format("%.4f		%.4f \n",
					-1 * solution.getObjective(0),
					-1 * solution.getObjective(0));
		}
		


			dke.finalValidation (dke.testingData, result.get(0));
	
			
			
	}
	
}

