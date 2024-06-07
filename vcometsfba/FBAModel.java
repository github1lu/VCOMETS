//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package edu.bu.segrelab.comets.fba;

import edu.bu.segrelab.comets.CometsConstants;
import edu.bu.segrelab.comets.Model;
import edu.bu.segrelab.comets.exception.ModelFileException;
import edu.bu.segrelab.comets.ui.DoubleField;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import org.apache.commons.lang3.ArrayUtils;

public class FBAModel extends Model implements CometsConstants {
    private int numRxns;
    private int numMetabs;
    private int numExch;
    private boolean runSuccess;
    private String modelID;
    private String ancestor;
    private String mutation;
    private int[] exch;
    private double[] exchDiffConsts;
    private String[] exchRxnNames;
    private String[] exchMetabNames;
    private String[] rxnNames;
    private String[] metabNames;
    private double[] baseLB;
    private double[] baseUB;
    private double[] baseExchLB;
    private double[] baseExchUB;
    private double[] exchKm;
    private double[] exchVmax;
    private double[] exchHillCoeff;
    private double[] exchAlpha;
    private double[] exchW;
    private double[][] lightAbsorption;
    private List<Signal> signals;
    private double flowDiffConst;
    private double growthDiffConst;
    private double elasticModulusConst;
    private double frictionConst;
    private double packedDensity;
    private double convectionDiffConst;
    private double convNonlinDiffZero;
    private double convNonlinDiffN;
    private double convNonlinDiffExponent;
    private double convNonlinDiffHillK;
    private double convNonlinDiffHillN;
    private double noiseVariance;
    private int[] objReactions;
    private double[] objWeights;
    private int biomassReaction;
    private int maintenanceReaction;
    private double maintenanceFlux;
    private int objStyle;
    private double defaultLB;
    private double defaultUB;
    private double defaultKm;
    private double defaultVmax;
    private double defaultHill;
    private double defaultAlpha;
    private double defaultW;
    private double defaultMetabDiffConst;
    private double genomeCost;
    private boolean active;
    public static final int MAXIMIZE_OBJECTIVE_FLUX = 0;
    public static final int MINIMIZE_OBJECTIVE_FLUX = 1;
    public static final int MAXIMIZE_TOTAL_FLUX = 2;
    public static final int MINIMIZE_TOTAL_FLUX = 3;
    public static final int MAX_OBJECTIVE_MIN_TOTAL = 4;
    public static final int MAX_OBJECTIVE_MAX_TOTAL = 5;
    public static final int MIN_OBJECTIVE_MIN_TOTAL = 6;
    public static final int MIN_OBJECTIVE_MAX_TOTAL = 7;
    private FBAOptimizer fbaOptimizer;
    public static final int GUROBI = 0;
    public static final int GLPK = 1;
    private ModelParametersPanel paramsPanel;
    private double neutralDriftSigma;
    private boolean neutralDrift;

    private FBAModel() {
        this.maintenanceReaction = 0;
        this.maintenanceFlux = 0.0;
        this.defaultLB = 0.0;
        this.defaultUB = 0.0;
        this.defaultKm = -1.0;
        this.defaultVmax = -1.0;
        this.defaultHill = -1.0;
        this.defaultAlpha = -1.0;
        this.defaultW = -1.0;
        this.defaultMetabDiffConst = 0.0;
        this.genomeCost = 0.0;
        this.neutralDriftSigma = 0.01;
        this.neutralDrift = false;
        this.runSuccess = false;
        this.objStyle = 0;
    }

    public FBAModel(double[][] m, double[] l, double[] u, int r, int optim) {
        this(m, l, u, new int[]{Math.abs(r)}, new double[]{r >= 0 ? 1.0 : -1.0}, Math.abs(r), optim);
    }

    public FBAModel(double[][] m, double[] l, double[] u, int[] objs, double[] objWt, int b, int optim) {
        this.maintenanceReaction = 0;
        this.maintenanceFlux = 0.0;
        this.defaultLB = 0.0;
        this.defaultUB = 0.0;
        this.defaultKm = -1.0;
        this.defaultVmax = -1.0;
        this.defaultHill = -1.0;
        this.defaultAlpha = -1.0;
        this.defaultW = -1.0;
        this.defaultMetabDiffConst = 0.0;
        this.genomeCost = 0.0;
        this.neutralDriftSigma = 0.01;
        this.neutralDrift = false;
        this.runSuccess = false;
        this.objStyle = 4;
        switch (optim) {
            case 0:
                this.fbaOptimizer = new FBAOptimizerGurobi(m, l, u, objs, objWt);
                break;
            case 1:
                this.fbaOptimizer = new FBAOptimizerGLPK(m, l, u, objs);
        }

        Double mtb = m[m.length - 1][0];
        this.numMetabs = mtb.intValue();
        this.numRxns = 1;

        for(int i = 0; i < m.length; ++i) {
            if (m[i][1] > (double)this.numRxns) {
                Double k = m[i][1];
                this.numRxns = k.intValue();
            }
        }

        this.setBaseBounds(l, u);
        this.setObjectiveReactions(objs);
        this.setObjectiveWeights(objWt);
        this.setBiomassReaction(b);
        this.signals = new ArrayList();
    }

    public FBAModel(double[][] m, double[] l, double[] u, int[] r, double[] objWt, int b, int[] exch, double[] exchDiffConsts, double[] exchKm, double[] exchVmax, double[] exchHillCoeff, double[] exchAlpha, double[] exchW, double[][] lightAbsorption, String[] metabNames, String[] rxnNames, int objStyle, int optim) {
        this(m, l, u, r, objWt, b, optim);
        if (exch == null) {
            throw new IllegalArgumentException("There must be an array of exchange reactions.");
        } else {
            this.exch = (int[])exch.clone();
            this.exchDiffConsts = exchDiffConsts;
            this.exchKm = exchKm;
            this.exchVmax = exchVmax;
            this.exchHillCoeff = exchHillCoeff;
            this.exchAlpha = exchAlpha;
            this.exchW = exchW;
            this.objStyle = objStyle;
            this.metabNames = metabNames;
            this.rxnNames = rxnNames;
            this.numExch = exch.length;
            this.exchRxnNames = new String[this.numExch];
            this.exchMetabNames = new String[this.numExch];
            int[] exch_tmp = (int[])exch.clone();

            int i;
            for(i = 0; i < m.length; ++i) {
                Double curr_rxn = m[i][1];
                Double curr_mtb = m[i][0];
                int is_exch = Arrays.binarySearch(exch_tmp, curr_rxn.intValue());
                if (is_exch >= 0) {
                    int rxn = curr_rxn.intValue();
                    int mtb = curr_mtb.intValue();
                    this.exchRxnNames[is_exch] = rxnNames[rxn - 1];
                    this.exchMetabNames[is_exch] = metabNames[mtb - 1];
                }
            }

            this.lightAbsorption = lightAbsorption;
            this.baseExchLB = new double[this.numExch];
            this.baseExchUB = new double[this.numExch];

            for(i = 0; i < this.numExch; ++i) {
                this.baseExchLB[i] = this.baseLB[exch[i] - 1];
                this.baseExchUB[i] = this.baseUB[exch[i] - 1];
            }

        }
    }

    public FBAModel(double[][] m, double[] l, double[] u, int[] r, double[] objWt, int[] exch, double[] exchDiffConsts, double[] exchKm, double[] exchVmax, double[] exchHillCoeff, double[] exchAlpha, double[] exchW, double[][] lightAbsorption, String[] metabNames, String[] rxnNames, int objStyle, int optim) {
        this(m, l, u, r, objWt, r[0], exch, exchDiffConsts, exchKm, exchVmax, exchHillCoeff, exchAlpha, exchW, lightAbsorption, metabNames, rxnNames, objStyle, optim);
    }

    public double getDefaultLB() {
        return this.defaultLB;
    }

    public double getDefaultUB() {
        return this.defaultUB;
    }

    public double getDefaultKm() {
        return this.defaultKm;
    }

    public double getDefaultVmax() {
        return this.defaultVmax;
    }

    public double getDefaultHill() {
        return this.defaultHill;
    }

    public double getDefaultAlpha() {
        return this.defaultAlpha;
    }

    public double getDefaultW() {
        return this.defaultW;
    }

    public double getDefaultMetabDiffConst() {
        return this.defaultMetabDiffConst;
    }

    public void setDefaultLB(double defLB) {
        this.defaultLB = defLB;
    }

    public void setDefaultUB(double defUB) {
        this.defaultUB = defUB;
    }

    public void setDefaultKm(double defKm) {
        this.defaultKm = defKm;
    }

    public void setDefaultVmax(double defVmax) {
        this.defaultVmax = defVmax;
    }

    public void setDefaultHill(double defHill) {
        this.defaultHill = defHill;
    }

    public void setDefaultAlpha(double defAlpha) {
        this.defaultAlpha = defAlpha;
    }

    public void setDefaultW(double defW) {
        this.defaultW = defW;
    }

    public void setDefaultMetabDiffConst(double defDiff) {
        for(int i = 0; i < this.exchDiffConsts.length; ++i) {
            if (this.exchDiffConsts[i] == this.defaultMetabDiffConst) {
                this.exchDiffConsts[i] = defDiff;
            }
        }

        this.defaultMetabDiffConst = defDiff;
    }

    public double[] getExchangeKm() {
        return this.exchKm;
    }

    public double[] getExchangeKmWithDefaults() {
        double[] res = new double[this.exchKm.length];

        for(int i = 0; i < this.exchKm.length; ++i) {
            res[i] = this.getExchangeKm(i);
        }

        return res;
    }

    public double getExchangeKm(int i) {
        double d = this.exchKm[i];
        if (d < 0.0) {
            d = this.defaultKm;
            if (d < 0.0) {
                d = FBAParameters.getDefaultKm();
            }
        }

        return d;
    }

    public void setExchangeKm(double[] exchKm) {
        if (this.numExch == exchKm.length) {
            this.exchKm = exchKm;
        }

    }

    public double[] getExchangeVmax() {
        return this.exchVmax;
    }

    public double[] getExchangeVmaxWithDefaults() {
        double[] res = new double[this.exchVmax.length];

        for(int i = 0; i < this.exchVmax.length; ++i) {
            res[i] = this.getExchangeVmax(i);
        }

        return res;
    }

    public double getExchangeVmax(int i) {
        double d = this.exchVmax[i];
        if (d < 0.0) {
            d = this.defaultVmax;
            if (d < 0.0) {
                d = FBAParameters.getDefaultVmax();
            }
        }

        return d;
    }

    public void setExchangeVmax(double[] exchVmax) {
        if (this.numExch == exchVmax.length) {
            this.exchVmax = exchVmax;
        }

    }

    public double[] getExchangeHillCoefficients() {
        return this.exchHillCoeff;
    }

    public double[] getExchangeHillCoefficientsWithDefaults() {
        double[] res = new double[this.exchHillCoeff.length];

        for(int i = 0; i < this.exchHillCoeff.length; ++i) {
            res[i] = this.getExchangeHillCoefficient(i);
        }

        return res;
    }

    public double getExchangeHillCoefficient(int i) {
        double d = this.exchHillCoeff[i];
        if (d < 0.0) {
            d = this.defaultHill;
            if (d < 0.0) {
                d = FBAParameters.getDefaultHill();
            }
        }

        return d;
    }

    public void setExchangeHillCoefficients(double[] exchHillCoeff) {
        if (this.numExch == exchHillCoeff.length) {
            this.exchHillCoeff = exchHillCoeff;
        }

    }

    public double[] getExchangeAlphaCoefficients() {
        return this.exchAlpha;
    }

    public double[] getExchangeAlphaCoefficientsWithDefaults() {
        double[] res = new double[this.exchAlpha.length];

        for(int i = 0; i < this.exchAlpha.length; ++i) {
            res[i] = this.getExchangeAlphaCoefficient(i);
        }

        return res;
    }

    public double getExchangeAlphaCoefficient(int i) {
        double d = this.exchAlpha[i];
        if (d < 0.0) {
            d = this.defaultAlpha;
            if (d < 0.0) {
                d = FBAParameters.getDefaultAlpha();
            }
        }

        return d;
    }

    public void setExchangeAlphaCoefficients(double[] exchAlphaCoeff) {
        if (this.numExch == exchAlphaCoeff.length) {
            this.exchAlpha = exchAlphaCoeff;
        }

    }

    public double[] getExchangeWCoefficients() {
        return this.exchW;
    }

    public double[] getExchangeWCoefficientsWithDefaults() {
        double[] res = new double[this.exchW.length];

        for(int i = 0; i < this.exchW.length; ++i) {
            res[i] = this.getExchangeWCoefficient(i);
        }

        return res;
    }

    public double getExchangeWCoefficient(int i) {
        double d = this.exchW[i];
        if (d < 0.0) {
            d = this.defaultW;
            if (d < 0.0) {
                d = FBAParameters.getDefaultW();
            }
        }

        return d;
    }

    public void setExchangeWCoefficients(double[] exchW) {
        if (this.numExch == exchW.length) {
            this.exchW = exchW;
        }

    }

    public double[] getExchangeMetaboliteDiffusionConstants() {
        return (double[])this.exchDiffConsts.clone();
    }

    public void setExchangeMetaboliteDiffusionConstants(double[] metabDiffConsts) {
        if (this.numExch == metabDiffConsts.length) {
            this.exchDiffConsts = metabDiffConsts;
        }

    }

    public double[][] getLightAbsorption() {
        return this.lightAbsorption;
    }

    public double[] getLightAbsorption(int i) {
        return this.lightAbsorption[i];
    }

    public void setLightAbsorption(double[][] lightAbsorption) {
        if (this.numExch == lightAbsorption.length) {
            this.lightAbsorption = lightAbsorption;
        }

    }

    public String getModelName() {
        return null;
    }

    public int setBaseBounds(double[] lb, double[] ub) {
        if (lb.length != ub.length) {
            return 0;
        } else {
            int i;
            for(i = 0; i < lb.length; ++i) {
                if (lb[i] > ub[i]) {
                    return 0;
                }
            }

            this.setBaseLowerBounds(lb);
            i = this.setBaseUpperBounds(ub);
            return i;
        }
    }

    public int setBaseLowerBounds(double[] lb) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            if (lb.length != this.numRxns) {
                return 0;
            } else {
                this.baseLB = lb;
                return 1;
            }
        } else {
            return 2;
        }
    }

    public int setBaseUpperBounds(double[] ub) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            if (this.numRxns != ub.length) {
                return 0;
            } else {
                this.baseUB = ub;
                return 1;
            }
        } else {
            return 2;
        }
    }

    public int setBaseExchLowerBounds(double[] lb) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            if (this.numExch != lb.length) {
                return 0;
            } else {
                this.baseExchLB = lb;
                return 1;
            }
        } else {
            return 2;
        }
    }

    public int setBaseExchUpperBounds(double[] ub) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            if (this.numExch != ub.length) {
                return 0;
            } else {
                this.baseExchUB = ub;
                return 1;
            }
        } else {
            return 2;
        }
    }

    public double[] getBaseLowerBounds() {
        return (double[])this.baseLB.clone();
    }

    public double[] getBaseUpperBounds() {
        return (double[])this.baseUB.clone();
    }

    public synchronized double[] getBaseExchLowerBounds() {
        return (double[])this.baseExchLB.clone();
    }

    public double[] getBaseExchUpperBounds() {
        return (double[])this.baseExchUB.clone();
    }

    public String[] getReactionNames() {
        return (String[])this.rxnNames.clone();
    }

    public void setReactionNames(String[] rxnNames) {
        if (this.numRxns == rxnNames.length) {
            this.rxnNames = rxnNames;
        }

    }

    public String[] getMetaboliteNames() {
        return (String[])this.metabNames.clone();
    }

    public void setMetaboliteNames(String[] metabNames) {
        if (this.numMetabs == metabNames.length) {
            this.metabNames = metabNames;
        }

    }

    public int setExchLowerBounds(double[] lb) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            if (lb.length != this.numExch) {
                return 0;
            } else {
                this.fbaOptimizer.setExchLowerBounds(this.exch, lb);
                return 1;
            }
        } else {
            return 2;
        }
    }

    public int setExchUpperBounds(double[] ub) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            if (ub.length != this.numExch) {
                return 0;
            } else {
                this.fbaOptimizer.setExchUpperBounds(this.exch, ub);
                return 1;
            }
        } else {
            return 2;
        }
    }

    public double[] getLowerBounds() {
        double[] l = new double[this.numRxns];
        l = this.fbaOptimizer.getLowerBounds(this.numRxns);
        return l;
    }

    public int setUpperBounds(double[] ub) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            this.fbaOptimizer.setUpperBounds(this.numRxns, ub);
            return 1;
        } else {
            return 2;
        }
    }

    public int setLowerBounds(double[] lb) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            this.fbaOptimizer.setLowerBounds(this.numRxns, lb);
            return 1;
        } else {
            return 2;
        }
    }

    public double[] getUpperBounds() {
        double[] u = new double[this.numRxns];
        u = this.fbaOptimizer.getUpperBounds(this.numRxns);
        return u;
    }

    public int setObjectiveStyle(int obj) {
        if (obj != 0 && obj != 1 && obj != 2 && obj != 3 && obj != 4 && obj != 5 && obj != 6 && obj != 7) {
            return 0;
        } else {
            this.objStyle = obj;
            return 1;
        }
    }

    public int getObjectiveStyle() {
        return this.objStyle;
    }

    public int[] getObjectiveIndexes() {
        return this.objReactions;
    }

    public int getObjectiveIndex() {
        return this.objReactions[0];
    }

    public int setObjectiveReaction(int r) {
        if (r >= 1 && r <= this.numRxns) {
            if (this.numMetabs != 0 && this.numRxns != 0) {
                this.fbaOptimizer.setObjectiveReaction(this.numRxns, r);
                this.objReactions = new int[]{r};
                return 1;
            } else {
                return 2;
            }
        } else {
            return 0;
        }
    }

    public int setObjectiveReactions(int[] objs) {
        this.objReactions = objs;
        return this.fbaOptimizer.setObjectiveReaction(this.numRxns, objs);
    }

    public int setObjectiveWeights(double[] objWt) {
        this.objWeights = objWt;
        return this.fbaOptimizer.setObjectiveWeights(objWt);
    }

    public List<Signal> getSignals() {
        return this.signals;
    }

    public void setSignals(List<Signal> signals) {
        this.signals = signals;
    }

    public synchronized int run() {
        this.runSuccess = false;
        if (this.numMetabs != 0 && this.numRxns != 0) {
            int ret = true;
            int ret = this.fbaOptimizer.run(this.objStyle);
            if (ret == 5) {
                this.runSuccess = true;
            }

            return ret;
        } else {
            return 2;
        }
    }

    public String getModelID() {
        return this.modelID;
    }

    public void setModelID(String model_id) {
        this.modelID = model_id;
    }

    public String getAncestor() {
        return this.ancestor;
    }

    public void setAncestor(String ancestor_id) {
        this.ancestor = ancestor_id;
    }

    public void setMutation(String mutation) {
        this.mutation = mutation;
    }

    public String getMutation() {
        return this.mutation;
    }

    public double[] getFluxes() {
        double[] v = new double[this.numRxns];
        if (this.runSuccess) {
            v = this.fbaOptimizer.getFluxes();
        }

        return v;
    }

    public double[] getExchangeFluxes() {
        double[] v = new double[this.numExch];
        if (this.runSuccess) {
            v = this.fbaOptimizer.getExchangeFluxes(this.exch);
        }

        return v;
    }

    public void setExchangeReactionNames(String[] rxnNames) {
        this.exchRxnNames = rxnNames;
    }

    public String[] getExchangeReactionNames() {
        return (String[])this.exchRxnNames.clone();
    }

    public void setExchangeMetaboliteNames(String[] metabNames) {
        this.exchMetabNames = metabNames;
    }

    public String[] getExchangeMetaboliteNames() {
        return (String[])this.exchMetabNames.clone();
    }

    public int[] getExchangeIndices() {
        return (int[])this.exch.clone();
    }

    public void setExchangeIndices(int[] exch) {
        if (this.numExch == exch.length) {
            this.exch = exch;
        }

    }

    public double[] getObjectiveSolutions() {
        return this.fbaOptimizer.getObjectiveSolutions(this.objReactions);
    }

    public double[] getObjectiveFluxSolution() {
        return this.fbaOptimizer.getObjectiveSolutions(this.objReactions);
    }

    public double getBiomassFluxSolution() {
        return this.fbaOptimizer.getObjectiveSolution(this.biomassReaction);
    }

    public double getMaintenanceFluxSolution() {
        return this.fbaOptimizer.getObjectiveSolution(this.maintenanceReaction);
    }

    public double getMinMaintenanceFlux() {
        return this.maintenanceFlux;
    }

    public boolean hasMaintenance() {
        return this.maintenanceReaction != 0;
    }

    public void setMaintenanceReaction(int rxn, double flux) {
        this.maintenanceReaction = rxn;
        this.maintenanceFlux = flux;
    }

    public int getFBAstatus() {
        return this.fbaOptimizer.getFBAstatus();
    }

    public static FBAModel loadModelFromFile(String filename) throws ModelFileException {
        try {
            int lineNum = 0;
            BufferedReader reader = new BufferedReader(new FileReader(filename));
            int lines_sparse_s = 0;
            BufferedReader reader_2 = new BufferedReader(new FileReader(filename));
            int numMets = 0;
            int numRxns = 0;
            double[][] S = null;
            double[] lb = null;
            double[] ub = null;
            int[] objs = new int[]{0};
            double[] objWt = new double[]{1.0};
            int bio = 0;
            int objSt = 0;
            int optim = 0;
            String[] metNames = null;
            String[] rxnNames = null;
            int[] exchRxns = null;
            int numExch = 0;
            double[] diffConsts = null;
            double[] exchAlpha = null;
            double[] exchW = null;
            double[] exchKm = null;
            double[] exchVmax = null;
            double[] exchHillCoeff = null;
            double[][] lightAbsorption = null;
            int maintRxn = 0;
            double maintFlux = 0.0;
            List<Signal> signals = new ArrayList();
            double defaultAlpha = -1.0;
            double defaultW = -1.0;
            double defaultKm = -1.0;
            double defaultVmax = -1.0;
            double defaultHill = -1.0;
            double defaultLB = -1000.0;
            double defaultUB = 1000.0;
            double defaultDiff = 1.0E-6;
            double elasticModulusConst = 1.0;
            double frictionConst = 1.0;
            double convDiffConst = 1.0;
            double convNonlinDiffZero = 1.0;
            double convNonlinDiffN = 1.0;
            double convNonlinDiffExponent = 1.0;
            double convNonlinDiffHillN = 10.0;
            double convNonlinDiffHillK = 0.9;
            double packDensity = 1.0;
            double noiseVariance = 0.0;
            double neutralDriftSigma = 0.0;
            boolean blockOpen = false;
            boolean neutralDrift = false;
            String line_2 = null;

            label1039:
            while((line_2 = reader_2.readLine()) != null) {
                if (line_2.contains("SMATRIX")) {
                    while(true) {
                        if ((line_2 = reader_2.readLine()) == null) {
                            break label1039;
                        }

                        ++lines_sparse_s;
                        if (line_2.contains("//")) {
                            break label1039;
                        }
                    }
                }
            }

            --lines_sparse_s;
            reader_2.close();
            String line = null;

            int i;
            label1028:
            while((line = reader.readLine()) != null) {
                line = line.trim();
                ++lineNum;
                String[] tokens = line.split("\\s+");
                String[] parsed;
                int j;
                int i;
                double u;
                if (tokens[0].equalsIgnoreCase("SMATRIX")) {
                    i = 0;
                    if (tokens.length != 3) {
                        reader.close();
                        throw new ModelFileException("The SMATRIX line should include the number of rows and columns of the Stoichiometric matrix on line " + lineNum);
                    }

                    numMets = Integer.parseInt(tokens[1]);
                    if (numMets <= 0) {
                        reader.close();
                        throw new ModelFileException("There must be at least one row (metabolite) in the Stoichiometric matrix.");
                    }

                    numRxns = Integer.parseInt(tokens[2]);
                    if (numRxns <= 0) {
                        reader.close();
                        throw new ModelFileException("There must be at least one column (reaction) in the Stoichiometric matrix.");
                    }

                    S = new double[lines_sparse_s][3];
                    parsed = null;
                    blockOpen = true;

                    while(true) {
                        String matLine;
                        do {
                            if ((matLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                ++lineNum;
                                blockOpen = false;
                                continue label1028;
                            }
                        } while(matLine.length() == 0);

                        String[] parsed = matLine.split("\\s+");
                        if (parsed.length != 3) {
                            reader.close();
                            throw new ModelFileException("Each line of the SMATRIX block should contain three elements - a row, column, and stoichiometric value for that element on line " + lineNum);
                        }

                        j = Integer.parseInt(parsed[0]);
                        if (j >= 1 && j <= numMets) {
                            i = Integer.parseInt(parsed[1]);
                            if (j >= 1 && j <= numMets) {
                                u = Double.parseDouble(parsed[2]);
                                S[i][0] = (double)j;
                                S[i][1] = (double)i;
                                S[i][2] = u;
                                ++i;
                                ++lineNum;
                                continue;
                            }

                            reader.close();
                            throw new ModelFileException("The second element of the SMATRIX block at line " + lineNum + " corresponds to the column, and should be between 1 and the number of columns specified.");
                        }

                        reader.close();
                        throw new ModelFileException("The first element of the SMATRIX block at line " + lineNum + " corresponds to the row, and should be between 1 and the number of rows specified.");
                    }
                } else {
                    String[] parsed;
                    int rxn;
                    String hillLine;
                    double hill;
                    if (tokens[0].equalsIgnoreCase("BOUNDS")) {
                        if (numRxns <= 0) {
                            reader.close();
                            throw new ModelFileException("The stoichiometric matrix should be loaded before the upper and lower bounds at line " + lineNum);
                        }

                        if (tokens.length != 3) {
                            reader.close();
                            throw new ModelFileException("The BOUNDS line should contain default lower and upper bound values, in that order at line " + lineNum);
                        }

                        defaultLB = Double.parseDouble(tokens[1]);
                        defaultUB = Double.parseDouble(tokens[2]);
                        if (defaultLB > defaultUB) {
                            reader.close();
                            throw new ModelFileException("The default lower bound should be LESS than the default upper bound at line " + lineNum);
                        }

                        lb = new double[numRxns];
                        ub = new double[numRxns];

                        for(i = 0; i < numRxns; ++i) {
                            lb[i] = defaultLB;
                            ub[i] = defaultUB;
                        }

                        parsed = null;
                        blockOpen = true;

                        while(true) {
                            do {
                                if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                    ++lineNum;
                                    blockOpen = false;
                                    continue label1028;
                                }

                                ++lineNum;
                            } while(hillLine.length() == 0);

                            parsed = hillLine.split("\\s+");
                            if (parsed.length != 3) {
                                reader.close();
                                throw new ModelFileException("There should be 3 elements on the BOUNDS line at file line " + lineNum + ": the reaction index (from 1 to N), the lower bound, and the upper bound.");
                            }

                            rxn = Integer.parseInt(parsed[0]);
                            if (rxn < 1 || rxn > numRxns) {
                                reader.close();
                                throw new ModelFileException("The reaction index in BOUNDS block line " + lineNum + " should be between 1 and " + numRxns);
                            }

                            hill = Double.parseDouble(parsed[1]);
                            u = Double.parseDouble(parsed[2]);
                            if (hill > u) {
                                reader.close();
                                throw new ModelFileException("The lower bound should be less than the upper bound on line " + lineNum);
                            }

                            lb[rxn - 1] = hill;
                            ub[rxn - 1] = u;
                        }
                    } else {
                        int rxn_num;
                        int p;
                        if (tokens[0].equalsIgnoreCase("OBJECTIVE")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the objective reaction at line " + lineNum);
                            }

                            parsed = null;
                            blockOpen = true;
                            rxn_num = Integer.parseInt(tokens[1]);
                            int[] rxns = new int[rxn_num];
                            double[] weights = new double[rxn_num];
                            i = 0;

                            while(!(hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                ++lineNum;
                                if (hillLine.length() != 0) {
                                    String[] parsed = hillLine.split("\\s+");
                                    if (parsed.length != 2) {
                                        reader.close();
                                        throw new ModelFileException("Each reaction must have a corresponding reaction weight.");
                                    }

                                    p = Integer.parseInt(parsed[0]);
                                    double w = Double.parseDouble(parsed[1]);
                                    rxns[i] = Math.abs(p);
                                    weights[i++] = p >= 0 ? w : -w;
                                }
                            }

                            objs = rxns;
                            objWt = weights;
                            ++lineNum;
                            blockOpen = false;
                        } else if (tokens[0].equalsIgnoreCase("OBJECTIVE_STYLE")) {
                            parsed = null;
                            blockOpen = true;

                            while(!(hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                ++lineNum;
                                if (hillLine.length() != 0) {
                                    parsed = hillLine.split("\\s+");
                                    if (parsed.length != 1) {
                                        reader.close();
                                        throw new ModelFileException("There should be just 1 element for the objective style line - the name of the objective style.");
                                    }

                                    if (parsed[0].equalsIgnoreCase("MAXIMIZE_OBJECTIVE_FLUX")) {
                                        objSt = 0;
                                    } else if (parsed[0].equalsIgnoreCase("MINIMIZE_OBJECTIVE_FLUX")) {
                                        objSt = 1;
                                    } else if (parsed[0].equalsIgnoreCase("MAXIMIZE_TOTAL_FLUX")) {
                                        objSt = 2;
                                    } else if (parsed[0].equalsIgnoreCase("MINIMIZE_TOTAL_FLUX")) {
                                        objSt = 3;
                                    } else if (parsed[0].equalsIgnoreCase("MAX_OBJECTIVE_MIN_TOTAL")) {
                                        objSt = 4;
                                    } else if (parsed[0].equalsIgnoreCase("MAX_OBJECTIVE_MAX_TOTAL")) {
                                        objSt = 5;
                                    } else if (parsed[0].equalsIgnoreCase("MIN_OBJECTIVE_MIN_TOTAL")) {
                                        objSt = 6;
                                    } else {
                                        if (!parsed[0].equalsIgnoreCase("MIN_OBJECTIVE_MAX_TOTAL")) {
                                            reader.close();
                                            throw new ModelFileException("Wrong OBJECTIVE_STYLE input value in model file.");
                                        }

                                        objSt = 7;
                                    }
                                }
                            }

                            ++lineNum;
                            blockOpen = false;
                        } else if (tokens[0].equalsIgnoreCase("BIOMASS")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the biomass reaction at line " + lineNum);
                            }

                            parsed = null;
                            blockOpen = true;

                            while(true) {
                                do {
                                    if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                        ++lineNum;
                                        blockOpen = false;
                                        continue label1028;
                                    }

                                    ++lineNum;
                                } while(hillLine.length() == 0);

                                parsed = hillLine.split("\\s+");
                                if (parsed.length != 1) {
                                    reader.close();
                                    throw new ModelFileException("There should be just 1 element for the biomass line - the index of the reaction.");
                                }

                                rxn = Integer.parseInt(parsed[0]);
                                if (rxn < 1 || rxn > numRxns) {
                                    reader.close();
                                    throw new ModelFileException("The reaction index in BIOMASS block line " + lineNum + " should be between 1 and " + numRxns);
                                }

                                bio = rxn;
                            }
                        } else if (tokens[0].equalsIgnoreCase("OPTIMIZER")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The OPTIMIZER should be followed only by the optimizer value. " + lineNum);
                            }

                            optim = -1;
                            if (tokens[1].equalsIgnoreCase("GUROBI")) {
                                optim = 0;
                            } else if (tokens[1].equalsIgnoreCase("GLPK")) {
                                optim = 1;
                            }
                        } else if (tokens[0].equalsIgnoreCase("neutralDrift")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The neutralDrift should be followed only by the value true or false at line " + lineNum);
                            }

                            neutralDrift = Boolean.parseBoolean(tokens[1]);
                        } else if (tokens[0].equalsIgnoreCase("neutralDriftSigma")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The neutralDriftSigma should be followed only by the value at line " + lineNum);
                            }

                            neutralDriftSigma = Double.parseDouble(tokens[1]);
                            if (neutralDriftSigma <= 0.0) {
                                reader.close();
                                throw new ModelFileException("The neutral drift sigma value given at line " + lineNum + "should be positive.");
                            }
                        } else if (tokens[0].equalsIgnoreCase("elasticModulus")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The ElasticModulus should be followed only by the modulus value at line " + lineNum);
                            }

                            elasticModulusConst = Double.parseDouble(tokens[1]);
                            if (elasticModulusConst < 0.0) {
                                reader.close();
                                throw new ModelFileException("The elastic modulus value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("packedDensity")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The packedDensity should be followed only by the value at line " + lineNum);
                            }

                            packDensity = Double.parseDouble(tokens[1]);
                            if (packDensity < 0.0) {
                                reader.close();
                                throw new ModelFileException("The packedDensity value given at line " + lineNum + "should be > 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("noiseVariance")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The noiseVariance should be followed only by the value at line " + lineNum);
                            }

                            noiseVariance = Double.parseDouble(tokens[1]);
                            if (noiseVariance < 0.0) {
                                reader.close();
                                throw new ModelFileException("The noiseVariance value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("frictionConstant")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The FrictionConstant should be followed only by its value at line " + lineNum);
                            }

                            frictionConst = Double.parseDouble(tokens[1]);
                            if (frictionConst <= 0.0) {
                                reader.close();
                                throw new ModelFileException("The frictionConstant value given at line " + lineNum + "should be > 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("convDiffConstant")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The convDiffConstant should be followed only by its value at line " + lineNum);
                            }

                            convDiffConst = Double.parseDouble(tokens[1]);
                            if (convDiffConst < 0.0) {
                                reader.close();
                                throw new ModelFileException("The convDiffConstant value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("convNonlinDiffZero")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffZero should be followed only by its value at line " + lineNum);
                            }

                            convNonlinDiffZero = Double.parseDouble(tokens[1]);
                            if (convDiffConst < 0.0) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffZero value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("convNonlinDiffN")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffN should be followed only by its value at line " + lineNum);
                            }

                            convNonlinDiffN = Double.parseDouble(tokens[1]);
                            if (convDiffConst < 0.0) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffN value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("convNonlinDiffExponent")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffExponent should be followed only by its value at line " + lineNum);
                            }

                            convNonlinDiffExponent = Double.parseDouble(tokens[1]);
                            if (convDiffConst < 0.0) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffExponent value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("convNonlinDiffHillK")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffHillK should be followed only by its value at line " + lineNum);
                            }

                            convNonlinDiffHillK = Double.parseDouble(tokens[1]);
                            if (convNonlinDiffHillK < 0.0) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffHillK value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("convNonlinDiffHillN")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffHillN should be followed only by its value at line " + lineNum);
                            }

                            convNonlinDiffHillN = Double.parseDouble(tokens[1]);
                            if (convNonlinDiffHillN < 0.0) {
                                reader.close();
                                throw new ModelFileException("The convNonlinDiffHillN value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("METABOLITE_NAMES")) {
                            if (numMets <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the metabolite names at line " + lineNum);
                            }

                            metNames = new String[numMets];
                            parsed = null;
                            rxn_num = 0;
                            blockOpen = true;

                            while(!(hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                ++lineNum;
                                if (hillLine.length() != 0) {
                                    if (rxn_num >= numMets) {
                                        reader.close();
                                        throw new ModelFileException("There must be one name for each metabolite, on each line of the METABOLITE_NAMES block. There's at least one extra at line " + lineNum);
                                    }

                                    metNames[rxn_num] = hillLine;
                                    ++rxn_num;
                                }
                            }

                            ++lineNum;
                            blockOpen = false;
                            if (rxn_num != numMets) {
                                reader.close();
                                throw new ModelFileException("There must be one name for each metabolite, on each line of the METABOLITE_NAMES block. There are apparently " + (numMets - rxn_num) + " names missing.");
                            }
                        } else if (tokens[0].equalsIgnoreCase("REACTION_NAMES")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the reaction names at line " + lineNum);
                            }

                            rxnNames = new String[numRxns];
                            parsed = null;
                            rxn_num = 0;
                            blockOpen = true;

                            while(!(hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                ++lineNum;
                                if (hillLine.length() != 0) {
                                    if (rxn_num >= numRxns) {
                                        reader.close();
                                        throw new ModelFileException("There must be one name for each reaction, on each line of the REACTION_NAMES block. There's at least one extra at line " + lineNum);
                                    }

                                    rxnNames[rxn_num] = hillLine;
                                    ++rxn_num;
                                }
                            }

                            ++lineNum;
                            blockOpen = false;
                            if (rxn_num != numRxns) {
                                reader.close();
                                throw new ModelFileException("There must be one name for each reaction, on each line of the REACTION_NAMES block. There are apparently " + (numRxns - rxn_num) + " names missing.");
                            }
                        } else if (tokens[0].equalsIgnoreCase("EXCHANGE_REACTIONS")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the exchange reaction list at line " + lineNum);
                            }

                            parsed = null;
                            blockOpen = true;

                            while(true) {
                                do {
                                    if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                        continue label1028;
                                    }

                                    ++lineNum;
                                } while(hillLine.length() == 0);

                                parsed = hillLine.split("\\s+");
                                if (parsed.length > numRxns) {
                                    reader.close();
                                    throw new ModelFileException("There should be, at most, " + numRxns + " values in the EXCHANGE_REACTIONS block. Looks like there's " + parsed.length + " instead");
                                }

                                if (parsed.length == 0) {
                                    exchRxns = new int[0];
                                } else {
                                    Set<Integer> exchSet = new HashSet();

                                    for(j = 0; j < parsed.length; ++j) {
                                        i = Integer.parseInt(parsed[j]);
                                        if (i < 1 || i > numRxns) {
                                            reader.close();
                                            throw new ModelFileException("Each exchange reaction should be between 1 and " + numRxns + " at line " + lineNum);
                                        }

                                        exchSet.add(i);
                                    }

                                    exchRxns = new int[exchSet.size()];
                                    Iterator<Integer> it = exchSet.iterator();

                                    for(i = 0; it.hasNext(); ++i) {
                                        exchRxns[i] = (Integer)it.next();
                                    }

                                    Arrays.sort(exchRxns);
                                }

                                ++lineNum;
                                blockOpen = false;
                                numExch = exchRxns.length;
                            }
                        } else if (tokens[0].equalsIgnoreCase("MAINTENANCE")) {
                            blockOpen = true;
                            parsed = null;

                            while(!(hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                ++lineNum;
                                if (hillLine.length() != 0) {
                                    parsed = hillLine.split("\\s+");
                                    if (parsed.length != 2) {
                                        reader.close();
                                        throw new ModelFileException("Please enter a maintenance reaction and its corresponding flux");
                                    }

                                    maintRxn = Integer.parseInt(parsed[0]);
                                    maintFlux = Double.parseDouble(parsed[1]);
                                }
                            }
                        } else if (tokens[0].equalsIgnoreCase("ALPHA_VALUES")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the alpha values at line " + lineNum);
                            }

                            if (exchRxns == null) {
                                reader.close();
                                throw new ModelFileException("The list of exchange reactions should be loaded before the alpha values at line " + lineNum);
                            }

                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The ALPHA_VALUES block header should be followed only by the default alpha value at line " + lineNum);
                            }

                            defaultAlpha = Double.parseDouble(tokens[1]);
                            if (defaultAlpha <= 0.0) {
                                reader.close();
                                throw new ModelFileException("The default alpha value given at line " + lineNum + "should be > 0");
                            }

                            exchAlpha = new double[numExch];

                            for(i = 0; i < numExch; ++i) {
                                exchAlpha[i] = -1.0;
                            }

                            parsed = null;
                            blockOpen = true;

                            while(true) {
                                do {
                                    if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                        ++lineNum;
                                        blockOpen = false;
                                        continue label1028;
                                    }

                                    ++lineNum;
                                } while(hillLine.length() == 0);

                                parsed = hillLine.split("\\s+");
                                if (parsed.length != 2) {
                                    reader.close();
                                    throw new ModelFileException("There should be 2 elements on each line of the ALPHA_VALUES block at line " + lineNum + ": the exchange reaction index (from 1 to " + numExch + ") and the alpha value of that reaction.");
                                }

                                rxn = Integer.parseInt(parsed[0]);
                                if (rxn < 1 || rxn > numExch) {
                                    reader.close();
                                    throw new ModelFileException("The reaction index in ALPHA_VALUES block line " + lineNum + " should be between 1 and " + numExch);
                                }

                                hill = Double.parseDouble(parsed[1]);
                                if (hill <= 0.0) {
                                    reader.close();
                                    throw new ModelFileException("The alpha value on line " + lineNum + " should be > 0");
                                }

                                exchAlpha[rxn - 1] = hill;
                            }
                        } else if (tokens[0].equalsIgnoreCase("W_VALUES")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the W values at line " + lineNum);
                            }

                            if (exchRxns == null) {
                                reader.close();
                                throw new ModelFileException("The list of exchange reactions should be loaded before the W values at line " + lineNum);
                            }

                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The W_VALUES block header should be followed only by the default W value at line " + lineNum);
                            }

                            defaultW = Double.parseDouble(tokens[1]);
                            if (defaultW <= 0.0) {
                                reader.close();
                                throw new ModelFileException("The default W value given at line " + lineNum + "should be > 0");
                            }

                            exchW = new double[numExch];

                            for(i = 0; i < numExch; ++i) {
                                exchW[i] = -1.0;
                            }

                            parsed = null;
                            blockOpen = true;

                            while(true) {
                                do {
                                    if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                        ++lineNum;
                                        blockOpen = false;
                                        continue label1028;
                                    }

                                    ++lineNum;
                                } while(hillLine.length() == 0);

                                parsed = hillLine.split("\\s+");
                                if (parsed.length != 2) {
                                    reader.close();
                                    throw new ModelFileException("There should be 2 elements on each line of the W_VALUES block at line " + lineNum + ": the exchange reaction index (from 1 to " + numExch + ") and the W value of that reaction.");
                                }

                                rxn = Integer.parseInt(parsed[0]);
                                if (rxn < 1 || rxn > numExch) {
                                    reader.close();
                                    throw new ModelFileException("The reaction index in W_VALUES block line " + lineNum + " should be between 1 and " + numExch);
                                }

                                hill = Double.parseDouble(parsed[1]);
                                if (hill <= 0.0) {
                                    reader.close();
                                    throw new ModelFileException("The W value on line " + lineNum + " should be > 0");
                                }

                                exchW[rxn - 1] = hill;
                            }
                        } else if (tokens[0].equalsIgnoreCase("KM_VALUES")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the Km values at line " + lineNum);
                            }

                            if (exchRxns == null) {
                                reader.close();
                                throw new ModelFileException("The list of exchange reactions should be loaded before the Km values at line " + lineNum);
                            }

                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The KM_VALUES block header should be followed only by the default Km value at line " + lineNum);
                            }

                            defaultKm = Double.parseDouble(tokens[1]);
                            if (defaultKm <= 0.0) {
                                reader.close();
                                throw new ModelFileException("The default Km value given at line " + lineNum + "should be > 0");
                            }

                            exchKm = new double[numExch];

                            for(i = 0; i < numExch; ++i) {
                                exchKm[i] = -1.0;
                            }

                            parsed = null;
                            blockOpen = true;

                            while(true) {
                                do {
                                    if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                        ++lineNum;
                                        blockOpen = false;
                                        continue label1028;
                                    }

                                    ++lineNum;
                                } while(hillLine.length() == 0);

                                parsed = hillLine.split("\\s+");
                                if (parsed.length != 2) {
                                    reader.close();
                                    throw new ModelFileException("There should be 2 elements on each line of the KM_VALUES block at line " + lineNum + ": the exchange reaction index (from 1 to " + numExch + ") and the Km value of that reaction.");
                                }

                                rxn = Integer.parseInt(parsed[0]);
                                if (rxn < 1 || rxn > numExch) {
                                    reader.close();
                                    throw new ModelFileException("The reaction index in KM_VALUES block line " + lineNum + " should be between 1 and " + numExch);
                                }

                                hill = Double.parseDouble(parsed[1]);
                                if (hill <= 0.0) {
                                    reader.close();
                                    throw new ModelFileException("The Km value on line " + lineNum + " should be > 0");
                                }

                                exchKm[rxn - 1] = hill;
                            }
                        } else if (tokens[0].equalsIgnoreCase("VMAX_VALUES")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the Vmax values at line " + lineNum);
                            }

                            if (exchRxns == null) {
                                reader.close();
                                throw new ModelFileException("The list of exchange reactions should be loaded before the Vmax values at line " + lineNum);
                            }

                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The VMAX_VALUES block header should be followed only by the default Vmax value at line " + lineNum);
                            }

                            defaultVmax = Double.parseDouble(tokens[1]);
                            if (defaultVmax <= 0.0) {
                                reader.close();
                                throw new ModelFileException("The default Vmax value given at line " + lineNum + "should be > 0");
                            }

                            exchVmax = new double[numExch];

                            for(i = 0; i < numExch; ++i) {
                                exchVmax[i] = -1.0;
                            }

                            parsed = null;
                            blockOpen = true;

                            while(true) {
                                do {
                                    if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                        ++lineNum;
                                        blockOpen = false;
                                        continue label1028;
                                    }

                                    ++lineNum;
                                    parsed = hillLine.split("\\s+");
                                } while(hillLine.length() == 0);

                                if (parsed.length != 2) {
                                    reader.close();
                                    throw new ModelFileException("There should be 2 elements on each line of the VMAX_VALUES block at line " + lineNum + ": the exchange reaction index (from 1 to " + numExch + ") and the Vmax value of that reaction.");
                                }

                                rxn = Integer.parseInt(parsed[0]);
                                if (rxn < 1 || rxn > numExch) {
                                    reader.close();
                                    throw new ModelFileException("The reaction index in VMAX_VALUES block line " + lineNum + " should be between 1 and " + numExch);
                                }

                                hill = Double.parseDouble(parsed[1]);
                                if (hill <= 0.0) {
                                    reader.close();
                                    throw new ModelFileException("The vMax value on line " + lineNum + " should be > 0");
                                }

                                exchVmax[rxn - 1] = hill;
                            }
                        } else if (tokens[0].equalsIgnoreCase("HILL_COEFFICIENTS")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the Hill coefficients at line " + lineNum);
                            }

                            if (exchRxns == null) {
                                reader.close();
                                throw new ModelFileException("The list of exchange reactions should be loaded before the Hill coefficients at line " + lineNum);
                            }

                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The HILL_COEFFICIENTS block header should be followed only by the default Hill coefficient at line " + lineNum);
                            }

                            defaultHill = Double.parseDouble(tokens[1]);
                            if (defaultHill < 0.0) {
                                reader.close();
                                throw new ModelFileException("The default Hill coefficient given at line " + lineNum + "should be >= 0");
                            }

                            exchHillCoeff = new double[numExch];

                            for(i = 0; i < numExch; ++i) {
                                exchHillCoeff[i] = -1.0;
                            }

                            parsed = null;
                            blockOpen = true;

                            while(true) {
                                do {
                                    if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                        ++lineNum;
                                        blockOpen = false;
                                        continue label1028;
                                    }

                                    ++lineNum;
                                } while(hillLine.length() == 0);

                                parsed = hillLine.split("\\s+");
                                if (parsed.length != 2) {
                                    reader.close();
                                    throw new ModelFileException("There should be 2 elements on each line of the HILL_COEFFICIENTS block at line " + lineNum + ": the exchange reaction index (from 1 to " + numExch + ") and the Hill coefficient of that reaction.");
                                }

                                rxn = Integer.parseInt(parsed[0]);
                                if (rxn < 1 || rxn > numExch) {
                                    reader.close();
                                    throw new ModelFileException("The reaction index in HILL_COEFFICIENTS block line " + lineNum + " should be between 1 and " + numExch);
                                }

                                hill = Double.parseDouble(parsed[1]);
                                if (hill < 0.0) {
                                    reader.close();
                                    throw new ModelFileException("The Hill coefficient on line " + lineNum + " should be >= 0");
                                }

                                exchHillCoeff[rxn - 1] = hill;
                            }
                        } else if (tokens[0].equalsIgnoreCase("neutralDriftParameter")) {
                            if (tokens.length != 2) {
                                reader.close();
                                throw new ModelFileException("The neutralDriftParameter should be followed only by its value at line " + lineNum);
                            }

                            neutralDriftSigma = Double.parseDouble(tokens[1]);
                            if (neutralDriftSigma < 0.0) {
                                reader.close();
                                throw new ModelFileException("The neutralDriftSigma value given at line " + lineNum + "should be => 0");
                            }
                        } else if (tokens[0].equalsIgnoreCase("LIGHT")) {
                            if (numRxns <= 0) {
                                reader.close();
                                throw new ModelFileException("The stoichiometric matrix should be loaded before the Light coefficients at line " + lineNum);
                            }

                            if (exchRxns == null) {
                                reader.close();
                                throw new ModelFileException("The list of exchange reactions should be loaded before the Light coefficients at line " + lineNum);
                            }

                            lightAbsorption = new double[numExch][2];

                            for(i = 0; i < numExch; ++i) {
                                for(rxn_num = 0; rxn_num < 2; ++rxn_num) {
                                    lightAbsorption[i][rxn_num] = 0.0;
                                }
                            }

                            parsed = null;
                            blockOpen = true;

                            label978:
                            while(true) {
                                do {
                                    if ((hillLine = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                        ++lineNum;
                                        blockOpen = false;
                                        continue label1028;
                                    }

                                    ++lineNum;
                                    parsed = hillLine.split("\\s+");
                                } while(hillLine.length() == 0);

                                if (parsed.length != 3) {
                                    reader.close();
                                    throw new ModelFileException("There should be 3 elements on each line of the LIGHT block at line " + lineNum + ": the exchange reaction index (from 1 to " + numExch + ") and the absorption coefficient of that reaction.");
                                }

                                rxn = Integer.parseInt(parsed[0]);
                                if (rxn >= 1 && rxn <= numExch) {
                                    j = 0;

                                    while(true) {
                                        if (j >= 2) {
                                            continue label978;
                                        }

                                        lightAbsorption[rxn - 1][j] = Double.parseDouble(parsed[j + 1]);
                                        if (lightAbsorption[rxn - 1][j] < 0.0 || lightAbsorption[rxn - 1][j] > 1.0) {
                                            reader.close();
                                            throw new ModelFileException("The absorption value on line " + lineNum + " should be between 0 and 1");
                                        }

                                        ++j;
                                    }
                                }

                                reader.close();
                                throw new ModelFileException("The reaction index in LIGHT block line " + lineNum + " should be between 1 and " + numExch);
                            }
                        } else if (tokens[0].equalsIgnoreCase("MET_REACTION_SIGNAL")) {
                            while(!(line = reader.readLine().trim()).equalsIgnoreCase("//")) {
                                parsed = line.split("\\s+");
                                if (parsed.length < 5) {
                                    reader.close();
                                    throw new ModelFileException("There must be at least five values given for each MET_REACTION_SIGNAL:\nrxn exch bound A K B,\nline num: " + lineNum);
                                }

                                rxn_num = -1;
                                if (!parsed[0].toLowerCase().equals("death")) {
                                    rxn_num = Integer.parseInt(parsed[0]);
                                    if (rxn_num > numRxns) {
                                        reader.close();
                                        throw new ModelFileException("first argument in MET_REACTION_SIGNAL must be the number of a reaction, < # reactions in S matrix. line num: " + lineNum);
                                    }
                                }

                                rxn = Integer.parseInt(parsed[1]);
                                if (rxn > numExch) {
                                    reader.close();
                                    throw new ModelFileException("second argument in MET_REACTION_SIGNAL must be the number of an exchange metabolite, < # exchange mets in model list. line num: " + lineNum);
                                }

                                String bound = parsed[2];
                                if (!bound.equalsIgnoreCase("lb") && !bound.equalsIgnoreCase("ub") && !bound.equalsIgnoreCase("consume_met") && !bound.equalsIgnoreCase("met_unchanged")) {
                                    reader.close();
                                    throw new ModelFileException("third argument in MET_REACTION_SIGNAL must be the string ub, lb or consume_met,met_unchanged designating the affected bound (or whether the metabolite is consumed for death-causing toxins). line num" + lineNum);
                                }

                                String function = parsed[3].toLowerCase();
                                double[] parameters = new double[parsed.length - 4];

                                for(p = 4; p < parsed.length; ++p) {
                                    parameters[p - 4] = Double.parseDouble(parsed[p]);
                                }

                                if (rxn_num != -1) {
                                    if (bound.equalsIgnoreCase("lb")) {
                                        signals.add(new Signal(true, false, false, rxn_num, rxn, function, parameters));
                                    } else {
                                        signals.add(new Signal(false, true, false, rxn_num, rxn, function, parameters));
                                    }
                                } else if (bound.equalsIgnoreCase("consume_met")) {
                                    signals.add(new Signal(false, false, true, rxn_num, rxn, function, parameters));
                                } else {
                                    signals.add(new Signal(false, false, false, rxn_num, rxn, function, parameters));
                                }
                            }
                        }
                    }
                }
            }

            reader.close();
            if (blockOpen) {
                throw new ModelFileException("Each data block is expected to end with '//' on a single line.");
            } else if (S == null) {
                throw new ModelFileException("To make an FBA model, a Stoichiometric matrix MUST be initialized!");
            } else if (lb != null && ub != null) {
                if (exchRxns == null) {
                    throw new ModelFileException("To make an FBA model, a set of exchange reactions MUST be initialized!");
                } else {
                    int i;
                    if (diffConsts == null) {
                        diffConsts = new double[numExch];

                        for(i = 0; i < numExch; ++i) {
                            diffConsts[i] = defaultDiff;
                        }
                    }

                    if (exchKm == null) {
                        exchKm = new double[numExch];

                        for(i = 0; i < numExch; ++i) {
                            exchKm[i] = defaultKm;
                        }
                    }

                    if (exchVmax == null) {
                        exchVmax = new double[numExch];

                        for(i = 0; i < numExch; ++i) {
                            exchVmax[i] = defaultVmax;
                        }
                    }

                    if (exchHillCoeff == null) {
                        exchHillCoeff = new double[numExch];

                        for(i = 0; i < numExch; ++i) {
                            exchHillCoeff[i] = defaultHill;
                        }
                    }

                    if (exchAlpha == null) {
                        exchAlpha = new double[numExch];

                        for(i = 0; i < numExch; ++i) {
                            exchAlpha[i] = defaultAlpha;
                        }
                    }

                    if (exchW == null) {
                        exchW = new double[numExch];

                        for(i = 0; i < numExch; ++i) {
                            exchW[i] = defaultW;
                        }
                    }

                    if (lightAbsorption == null) {
                        lightAbsorption = new double[numExch][2];

                        for(i = 0; i < numExch; ++i) {
                            for(i = 0; i < 2; ++i) {
                                lightAbsorption[i][i] = 0.0;
                            }
                        }
                    }

                    if (bio == 0) {
                        bio = objs[0];
                    }

                    FBAModel model = new FBAModel(S, lb, ub, objs, objWt, bio, exchRxns, diffConsts, exchKm, exchVmax, exchHillCoeff, exchAlpha, exchW, lightAbsorption, metNames, rxnNames, objSt, optim);
                    model.setDefaultAlpha(defaultAlpha);
                    model.setDefaultW(defaultW);
                    model.setDefaultHill(defaultHill);
                    model.setDefaultKm(defaultKm);
                    model.setDefaultVmax(defaultVmax);
                    model.setDefaultLB(defaultLB);
                    model.setDefaultUB(defaultUB);
                    model.setDefaultMetabDiffConst(defaultDiff);
                    model.setElasticModulusConstant(elasticModulusConst);
                    model.setFrictionConstant(frictionConst);
                    model.setConvDiffConstant(convDiffConst);
                    model.setConvNonlinDiffZero(convNonlinDiffZero);
                    model.setConvNonlinDiffN(convNonlinDiffN);
                    model.setConvNonlinDiffHillK(convNonlinDiffHillK);
                    model.setConvNonlinDiffHillN(convNonlinDiffHillN);
                    model.setConvNonlinDiffExponent(convNonlinDiffExponent);
                    model.setPackedDensity(packDensity);
                    model.setNoiseVariance(noiseVariance);
                    model.setSignals(signals);
                    model.setNeutralDrift(neutralDrift);
                    model.setNeutralDriftSigma(neutralDriftSigma);
                    model.setMaintenanceReaction(maintRxn, maintFlux);
                    model.setFileName(filename);
                    return model;
                }
            } else {
                throw new ModelFileException("To make an FBA model, a set of lower and upper bounds MUST be initialized!");
            }
        } catch (FileNotFoundException var82) {
            throw new ModelFileException(ModelFileException.FILE_NOT_FOUND, "Unable to find model file '" + filename + "': " + var82);
        } catch (IOException var83) {
            throw new ModelFileException(ModelFileException.IO_ERROR, "I/O error in model file '" + filename + "': " + var83);
        } catch (NumberFormatException var84) {
            throw new ModelFileException(ModelFileException.NUMBER_FORMAT_ERROR, "Number formatting error in model file '" + filename + "': " + var84);
        }
    }

    public String[] getMediaNames() {
        return this.getExchangeMetaboliteNames();
    }

    public double getGrowthDiffusionConstant() {
        return this.growthDiffConst;
    }

    public void setGrowthDiffusionConstant(double val) {
        if (val >= 0.0) {
            this.growthDiffConst = val;
        }

    }

    public double getFlowDiffusionConstant() {
        return this.flowDiffConst;
    }

    public void setFlowDiffusionConstant(double val) {
        if (val >= 0.0) {
            this.flowDiffConst = val;
        }

    }

    public double getPackedDensity() {
        return this.packedDensity;
    }

    public void setPackedDensity(double density) {
        this.packedDensity = density;
    }

    public double getElasticModulusConstant() {
        return this.elasticModulusConst;
    }

    public void setElasticModulusConstant(double val) {
        this.elasticModulusConst = val;
    }

    public double getConvDiffConstant() {
        return this.convectionDiffConst;
    }

    public void setConvDiffConstant(double val) {
        this.convectionDiffConst = val;
    }

    public double getConvNonlinDiffZero() {
        return this.convNonlinDiffZero;
    }

    public void setConvNonlinDiffZero(double val) {
        this.convNonlinDiffZero = val;
    }

    public double getConvNonlinDiffN() {
        return this.convNonlinDiffN;
    }

    public void setConvNonlinDiffN(double val) {
        this.convNonlinDiffN = val;
    }

    public double getConvNonlinDiffExponent() {
        return this.convNonlinDiffExponent;
    }

    public void setConvNonlinDiffExponent(double val) {
        this.convNonlinDiffExponent = val;
    }

    public double getConvNonlinDiffHillK() {
        return this.convNonlinDiffHillK;
    }

    public void setConvNonlinDiffHillK(double val) {
        this.convNonlinDiffHillK = val;
    }

    public double getConvNonlinDiffHillN() {
        return this.convNonlinDiffHillN;
    }

    public void setConvNonlinDiffHillN(double val) {
        this.convNonlinDiffHillN = val;
    }

    public double getNoiseVariance() {
        return this.noiseVariance;
    }

    public void setNoiseVariance(double variance) {
        this.noiseVariance = variance;
    }

    public double getFrictionConstant() {
        return this.frictionConst;
    }

    public void setFrictionConstant(double val) {
        this.frictionConst = val;
    }

    public double getNeutralDriftSigma() {
        return this.neutralDriftSigma;
    }

    public void setNeutralDriftSigma(double val) {
        this.neutralDriftSigma = val;
    }

    public boolean getNeutralDrift() {
        return this.neutralDrift;
    }

    public void setNeutralDrift(boolean val) {
        this.neutralDrift = val;
    }

    public int setObjectiveUpperBound(double ub) {
        int res = 1;
        int[] var4 = this.objReactions;
        int var5 = var4.length;

        for(int var6 = 0; var6 < var5; ++var6) {
            int objReaction = var4[var6];
            int t = this.fbaOptimizer.setObjectiveUpperBound(objReaction, ub);
            if (t != 1) {
                res = 0;
            }
        }

        return res;
    }

    public int setObjectiveUpperBounds(double[] ub) {
        int res = 1;

        for(int i = 0; i < ub.length; ++i) {
            int t = this.fbaOptimizer.setObjectiveUpperBound(this.objReactions[i], ub[i]);
            if (t != 1) {
                res = 0;
            }
        }

        return res;
    }

    public int setBiomassUpperBound(double ub) {
        return this.fbaOptimizer.setObjectiveUpperBound(this.biomassReaction, ub);
    }

    public int setObjectiveLowerBound(double lb) {
        int res = 1;
        int[] var4 = this.objReactions;
        int var5 = var4.length;

        for(int var6 = 0; var6 < var5; ++var6) {
            int objReaction = var4[var6];
            int t = this.fbaOptimizer.setObjectiveLowerBound(objReaction, lb);
            if (t != 1) {
                res = 0;
            }
        }

        return res;
    }

    public int setObjectiveLowerBounds(double[] lb) {
        int res = 1;

        for(int i = 0; i < lb.length; ++i) {
            int t = this.fbaOptimizer.setObjectiveLowerBound(this.objReactions[i], lb[i]);
            if (t != 1) {
                res = 0;
            }
        }

        return res;
    }

    public JComponent getInfoPanel() {
        JPanel panel = new JPanel();
        panel.add(new JLabel("Model info goes here."));
        return panel;
    }

    private void setNums(int numMetabs, int numRxns, int numExch) {
        this.numMetabs = numMetabs;
        this.numRxns = numRxns;
        this.numExch = numExch;
    }

    public FBAModel clone() {
        FBAModel modelCopy = new FBAModel();
        modelCopy.setNums(this.numMetabs, this.numRxns, this.numExch);
        modelCopy.fbaOptimizer = this.fbaOptimizer.clone();
        modelCopy.setBaseBounds(this.getBaseLowerBounds(), this.getBaseUpperBounds());
        modelCopy.setBaseExchLowerBounds(this.getBaseExchLowerBounds());
        modelCopy.setBaseExchUpperBounds(this.getBaseExchUpperBounds());
        modelCopy.setObjectiveReactions(this.getObjectiveIndexes());
        modelCopy.setBiomassReaction(this.getBiomassReaction());
        modelCopy.setExchangeIndices(this.getExchangeIndices());
        modelCopy.setExchangeKm(this.getExchangeKm());
        modelCopy.setExchangeVmax(this.getExchangeVmax());
        modelCopy.setExchangeHillCoefficients(this.getExchangeHillCoefficients());
        modelCopy.setExchangeAlphaCoefficients(this.getExchangeAlphaCoefficients());
        modelCopy.setExchangeWCoefficients(this.getExchangeWCoefficients());
        modelCopy.setExchangeMetaboliteDiffusionConstants(this.getExchangeMetaboliteDiffusionConstants());
        modelCopy.setMetaboliteNames(this.getMetaboliteNames());
        modelCopy.setReactionNames(this.getReactionNames());
        modelCopy.setExchangeMetaboliteNames(this.getExchangeMetaboliteNames());
        modelCopy.setExchangeReactionNames(this.getExchangeReactionNames());
        modelCopy.setDefaultAlpha(this.getDefaultAlpha());
        modelCopy.setDefaultW(this.getDefaultW());
        modelCopy.setDefaultHill(this.getDefaultHill());
        modelCopy.setDefaultKm(this.getDefaultKm());
        modelCopy.setDefaultVmax(this.getDefaultVmax());
        modelCopy.setDefaultLB(this.getDefaultLB());
        modelCopy.setDefaultUB(this.getDefaultUB());
        modelCopy.setDefaultMetabDiffConst(this.getDefaultMetabDiffConst());
        modelCopy.setActive(this.getActive());
        modelCopy.setObjectiveStyle(this.getObjectiveStyle());
        modelCopy.setFileName(this.getFileName());
        modelCopy.setElasticModulusConstant(this.getElasticModulusConstant());
        modelCopy.setFrictionConstant(this.getFrictionConstant());
        modelCopy.setConvDiffConstant(this.getConvDiffConstant());
        modelCopy.setPackedDensity(this.getPackedDensity());
        modelCopy.setNoiseVariance(this.getNoiseVariance());
        modelCopy.setLightAbsorption(this.getLightAbsorption());
        modelCopy.setSignals(this.getSignals());
        return modelCopy;
    }

    public JComponent getParametersPanel() {
        this.paramsPanel = new ModelParametersPanel(this);
        return this.paramsPanel;
    }

    public void applyParameters() {
        if (this.paramsPanel != null) {
            this.paramsPanel.updateModelParameters();
        }
    }

    public boolean activate(double activationRate) {
        if (!this.active) {
            Random random = new Random();
            double r = random.nextDouble();
            if (r < activationRate) {
                this.active = true;
            }
        }

        return this.active;
    }

    public boolean getActive() {
        return this.active;
    }

    public void setActive(boolean act) {
        this.active = act;
    }

    public int getBiomassReaction() {
        return this.biomassReaction;
    }

    public void setBiomassReaction(int biomassReaction) {
        this.biomassReaction = biomassReaction;
    }

    public int getTotalRxns() {
        double[] lBounds = this.getBaseLowerBounds();
        double[] uBounds = this.getBaseUpperBounds();
        int totalRxns = 0;

        for(int j = 0; j < lBounds.length; ++j) {
            if ((lBounds[j] != 0.0 || uBounds[j] != 0.0) && !ArrayUtils.contains(this.exch, j)) {
                ++totalRxns;
            }
        }

        return totalRxns;
    }

    public int getInactiveRxns() {
        double[] lBounds = this.getBaseLowerBounds();
        double[] uBounds = this.getBaseUpperBounds();
        int totalInactiveRxns = 0;

        for(int j = 0; j < lBounds.length; ++j) {
            if ((lBounds[j] == 0.0 || uBounds[j] == 0.0) && !ArrayUtils.contains(this.exch, j)) {
                ++totalInactiveRxns;
            }
        }

        return totalInactiveRxns;
    }

    public void mutateModel() {
        double[] lBounds = this.getBaseLowerBounds();
        double[] uBounds = this.getBaseUpperBounds();
        ArrayList<Integer> nonzeroRxns = new ArrayList();

        int j;
        for(j = 0; j < lBounds.length; ++j) {
            if ((lBounds[j] != 0.0 || uBounds[j] != 0.0) && !ArrayUtils.contains(this.exch, j)) {
                nonzeroRxns.add(j);
            }
        }

        j = (Integer)nonzeroRxns.get((new Random()).nextInt(nonzeroRxns.size()));
        this.setMutation("del_" + Integer.toString(j));
        lBounds[j] = 0.0;
        uBounds[j] = 0.0;
        this.setBaseLowerBounds(lBounds);
        this.setBaseUpperBounds(uBounds);
        this.fbaOptimizer.setLowerBounds(lBounds.length, lBounds);
        this.fbaOptimizer.setUpperBounds(uBounds.length, uBounds);
    }

    public void addReactionToModel() {
        double[] lBounds = this.getBaseLowerBounds();
        double[] uBounds = this.getBaseUpperBounds();
        ArrayList<Integer> nonzeroRxns = new ArrayList();

        int mutReaction;
        for(mutReaction = 0; mutReaction < lBounds.length; ++mutReaction) {
            if (lBounds[mutReaction] == 0.0 && uBounds[mutReaction] == 0.0 && !ArrayUtils.contains(this.exch, mutReaction)) {
                nonzeroRxns.add(mutReaction);
            }
        }

        if (nonzeroRxns.size() > 0) {
            mutReaction = (Integer)nonzeroRxns.get((new Random()).nextInt(nonzeroRxns.size()));
            this.setMutation("add_" + Integer.toString(mutReaction));
            uBounds[mutReaction] = 1000.0;
            this.setBaseUpperBounds(uBounds);
            this.fbaOptimizer.setLowerBounds(lBounds.length, lBounds);
            this.fbaOptimizer.setUpperBounds(uBounds.length, uBounds);
        }

    }

    public double getGenomeCost() {
        return this.genomeCost;
    }

    public void setGenomeCost(double ind_frac_cost) {
        double[] lBounds = this.getBaseLowerBounds();
        double[] uBounds = this.getBaseUpperBounds();
        int num_reactions = 0;

        for(int j = 0; j < lBounds.length; ++j) {
            if (lBounds[j] != 0.0 || uBounds[j] != 0.0) {
                ++num_reactions;
            }
        }

        this.genomeCost = (double)(num_reactions * num_reactions) * ind_frac_cost;
    }

    private class ModelParametersPanel extends JPanel {
        private static final long serialVersionUID = -6609420365194597487L;
        private int objRxnIndex;
        private int objStyle;
        private JComboBox rxnNamesBox;
        private ButtonGroup fbaObjGroup;
        private FBAModel model;
        private JRadioButton maxObjButton;
        private JRadioButton minObjButton;
        private JRadioButton maxFluxButton;
        private JRadioButton minFluxButton;
        private JRadioButton maxObjMinFluxButton;
        private JRadioButton minObjMinFluxButton;
        private JRadioButton maxObjMaxFluxButton;
        private JRadioButton minObjMaxFluxButton;
        private DoubleField flowConstField;
        private DoubleField growthConstField;
        private DoubleField elasticModulusField;
        private DoubleField frictionConstField;
        private DoubleField convDiffConstField;
        private DoubleField packedDensityField;
        private DoubleField noiseVarianceField;

        public ModelParametersPanel(FBAModel model) {
            this.model = model;
            this.objRxnIndex = model.getObjectiveIndex();
            String[] rxns = model.getReactionNames();
            this.rxnNamesBox = new JComboBox(rxns);
            this.rxnNamesBox.setSelectedIndex(this.objRxnIndex - 1);
            this.objStyle = model.getObjectiveStyle();
            this.fbaObjGroup = new ButtonGroup();
            this.maxObjButton = new JRadioButton("Maximize objective reaction");
            this.minObjButton = new JRadioButton("Minimize objective reaction");
            this.maxFluxButton = new JRadioButton("Maximize total flux");
            this.minFluxButton = new JRadioButton("Minimize total flux");
            this.maxObjMinFluxButton = new JRadioButton("Max objective / Min flux");
            this.minObjMinFluxButton = new JRadioButton("Min objective / Min flux");
            this.maxObjMaxFluxButton = new JRadioButton("Max objective / Max flux");
            this.minObjMaxFluxButton = new JRadioButton("Min objective / Max flux");
            this.maxObjMaxFluxButton.setEnabled(false);
            this.minObjMaxFluxButton.setEnabled(false);
            JLabel flowConstLabel = new JLabel("Flow diffusion constant (cm^2/s): ", 2);
            this.flowConstField = new DoubleField(model.getFlowDiffusionConstant(), 6, false);
            JLabel growthConstLabel = new JLabel("Growth diffusion constant (cm^2/s): ", 2);
            this.growthConstField = new DoubleField(model.getGrowthDiffusionConstant(), 6, false);
            JLabel elasticModulusLabel = new JLabel("Elastic modulus constant (Pa): ", 2);
            this.elasticModulusField = new DoubleField(model.getElasticModulusConstant(), 6, false);
            JLabel frictionConstLabel = new JLabel("Friction constant (cm^2/s): ", 2);
            this.frictionConstField = new DoubleField(model.getFrictionConstant(), 6, false);
            JLabel convDiffConstLabel = new JLabel("Conv. model diffusion constant (cm^2/s): ", 2);
            this.convDiffConstField = new DoubleField(model.getConvDiffConstant(), 6, false);
            JLabel packedDensityLabel = new JLabel("Packed Density (g/cm^2) or (g/cm^3): ", 2);
            this.packedDensityField = new DoubleField(model.getPackedDensity(), 6, false);
            JLabel noiseVarianceLabel = new JLabel("Noise variance: ", 2);
            this.noiseVarianceField = new DoubleField(model.getNoiseVariance(), 6, false);
            this.fbaObjGroup.add(this.maxObjButton);
            this.fbaObjGroup.add(this.minObjButton);
            this.fbaObjGroup.add(this.maxFluxButton);
            this.fbaObjGroup.add(this.minFluxButton);
            this.fbaObjGroup.add(this.maxObjMinFluxButton);
            this.fbaObjGroup.add(this.minObjMinFluxButton);
            this.fbaObjGroup.add(this.maxObjMaxFluxButton);
            this.fbaObjGroup.add(this.minObjMaxFluxButton);
            this.setSelectedObjectiveButton(this.objStyle);
            GridBagConstraints gbc = new GridBagConstraints();
            this.setLayout(new GridBagLayout());
            gbc.gridx = 0;
            gbc.gridy = 0;
            gbc.gridwidth = 2;
            gbc.anchor = 21;
            this.add(new JLabel("Objective Reaction"), gbc);
            ++gbc.gridy;
            this.add(this.rxnNamesBox, gbc);
            ++gbc.gridy;
            this.add(Box.createVerticalStrut(10), gbc);
            ++gbc.gridy;
            this.add(this.maxObjButton, gbc);
            ++gbc.gridy;
            this.add(this.minObjButton, gbc);
            ++gbc.gridy;
            this.add(this.maxFluxButton, gbc);
            ++gbc.gridy;
            this.add(this.minFluxButton, gbc);
            ++gbc.gridy;
            this.add(this.maxObjMinFluxButton, gbc);
            ++gbc.gridy;
            this.add(this.minObjMinFluxButton, gbc);
            ++gbc.gridy;
            ++gbc.gridy;
            gbc.gridwidth = 1;
            this.add(flowConstLabel, gbc);
            gbc.gridx = 1;
            this.add(this.flowConstField, gbc);
            ++gbc.gridy;
            gbc.gridx = 0;
            this.add(growthConstLabel, gbc);
            gbc.gridx = 1;
            this.add(this.growthConstField, gbc);
            ++gbc.gridy;
            gbc.gridx = 0;
            this.add(elasticModulusLabel, gbc);
            gbc.gridx = 1;
            this.add(this.elasticModulusField, gbc);
            ++gbc.gridy;
            gbc.gridx = 0;
            this.add(frictionConstLabel, gbc);
            gbc.gridx = 1;
            this.add(this.frictionConstField, gbc);
            ++gbc.gridy;
            gbc.gridx = 0;
            this.add(convDiffConstLabel, gbc);
            gbc.gridx = 1;
            this.add(this.convDiffConstField, gbc);
            ++gbc.gridy;
            gbc.gridx = 0;
            this.add(packedDensityLabel, gbc);
            gbc.gridx = 1;
            this.add(this.packedDensityField, gbc);
            ++gbc.gridy;
            gbc.gridx = 0;
            this.add(noiseVarianceLabel, gbc);
            gbc.gridx = 1;
            this.add(this.noiseVarianceField, gbc);
        }

        public void updateModelParameters() {
            this.model.setObjectiveReaction(this.rxnNamesBox.getSelectedIndex() + 1);
            this.model.setObjectiveStyle(this.getSelectedObjectiveStyle());
            this.model.setGrowthDiffusionConstant(this.growthConstField.getDoubleValue());
            this.model.setFlowDiffusionConstant(this.flowConstField.getDoubleValue());
            this.model.setElasticModulusConstant(this.elasticModulusField.getDoubleValue());
            this.model.setFrictionConstant(this.frictionConstField.getDoubleValue());
            this.model.setConvDiffConstant(this.convDiffConstField.getDoubleValue());
            this.model.setPackedDensity(this.packedDensityField.getDoubleValue());
            this.model.setNoiseVariance(this.noiseVarianceField.getDoubleValue());
        }

        private void setSelectedObjectiveButton(int objStyle) {
            this.fbaObjGroup.clearSelection();
            switch (objStyle) {
                case 0:
                    this.maxObjButton.setSelected(true);
                    break;
                case 1:
                    this.minObjButton.setSelected(true);
                    break;
                case 2:
                    this.maxFluxButton.setSelected(true);
                    break;
                case 3:
                    this.minFluxButton.setSelected(true);
                    break;
                case 4:
                    this.maxObjMinFluxButton.setSelected(true);
                    break;
                case 5:
                    this.maxObjMaxFluxButton.setSelected(true);
                    break;
                case 6:
                    this.minObjMinFluxButton.setSelected(true);
                    break;
                case 7:
                    this.minObjMaxFluxButton.setSelected(true);
            }

        }

        private int getSelectedObjectiveStyle() {
            if (this.maxObjButton.isSelected()) {
                return 0;
            } else if (this.minObjButton.isSelected()) {
                return 1;
            } else if (this.maxFluxButton.isSelected()) {
                return 2;
            } else if (this.minFluxButton.isSelected()) {
                return 3;
            } else if (this.maxObjMinFluxButton.isSelected()) {
                return 4;
            } else if (this.minObjMinFluxButton.isSelected()) {
                return 6;
            } else if (this.maxObjMaxFluxButton.isSelected()) {
                return 5;
            } else {
                return this.minObjMaxFluxButton.isSelected() ? 7 : 0;
            }
        }
    }

    private class ModelInfoPanel extends JPanel {
        private static final long serialVersionUID = -5975621566079865933L;

        private ModelInfoPanel() {
        }
    }
}
