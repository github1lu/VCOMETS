//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package edu.bu.segrelab.comets.fba;

import edu.bu.segrelab.comets.CometsConstants;

public abstract class FBAOptimizer implements CometsConstants {
    public FBAOptimizer() {
    }

    public abstract int run(int var1);

    public abstract int setExchLowerBounds(int[] var1, double[] var2);

    public abstract int setExchUpperBounds(int[] var1, double[] var2);

    public abstract int setLowerBounds(int var1, double[] var2);

    public abstract double[] getLowerBounds(int var1);

    public abstract int setUpperBounds(int var1, double[] var2);

    public abstract double[] getUpperBounds(int var1);

    public abstract int setObjectiveReaction(int var1, int var2);

    public abstract int setObjectiveReaction(int var1, int[] var2);

    public abstract double getObjectiveSolution(int var1);

    public abstract double getObjectiveFluxSolution(int var1);

    public abstract int getFBAstatus();

    public abstract int setObjectiveUpperBound(int var1, double var2);

    public abstract int setObjectiveLowerBound(int var1, double var2);

    public abstract double[] getFluxes();

    public abstract double[] getExchangeFluxes(int[] var1);

    public abstract FBAOptimizer clone();

    public abstract double[] getObjectiveSolutions(int[] var1);

    public abstract int setObjectiveWeights(double[] var1);
}
