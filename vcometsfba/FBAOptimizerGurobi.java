//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package edu.bu.segrelab.comets.fba;

import edu.bu.segrelab.comets.CometsConstants;
import gurobi.GRBConstr;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import gurobi.GRB.DoubleAttr;
import gurobi.GRB.IntAttr;
import gurobi.GRB.IntParam;
import gurobi.GRB.StringAttr;
import java.util.Arrays;
import java.util.Comparator;

public class FBAOptimizerGurobi extends FBAOptimizer implements CometsConstants {
    private GRBEnv env;
    private GRBModel model;
    private GRBVar[] rxnFluxes;
    private GRBLinExpr[] rxnExpressions;
    private double[] fluxesModel;
    private double[][] stoichMatrix;
    private GRBEnv envMin;
    private GRBModel modelMin;
    private GRBVar[] modelMinVars;
    private int nVars;
    private final String biomassConstraintName;
    private boolean runSuccess;
    private int numMetabs;
    private int numRxns;
    private int[] objReactions;
    private double[] objWeights;
    private String[] objConstraintNames;

    private FBAOptimizerGurobi() {
        this.biomassConstraintName = "biomassConstraint";
        this.runSuccess = false;

        try {
            this.env = new GRBEnv();
            this.env.set(IntParam.OutputFlag, 0);
            this.model = new GRBModel(this.env);
        } catch (GRBException var2) {
            System.out.println("Error in FBAOptimizerGurobi private constructor method");
            System.out.println("Error code: " + var2.getErrorCode() + ". " + var2.getMessage());
        }

    }

    public FBAOptimizerGurobi(double[][] m, double[] l, double[] u, int[] objIdxs, double[] objWeights) {
        this();
        this.sortByColumn(m, 0);
        this.initializeAbsModel(m, l, u, objIdxs, objWeights);
        this.stoichMatrix = m;
        this.fluxesModel = new double[l.length];
        double[] objective = new double[l.length];
        char[] types = new char[l.length];

        for(int i = 0; i < l.length; ++i) {
            objective[i] = 0.0;
            types[i] = 'C';
        }

        try {
            this.rxnFluxes = this.model.addVars(l, u, objective, types, (String[])null);
            this.model.update();
            Double mtb = m[m.length - 1][0];
            this.numMetabs = mtb.intValue();
            this.numRxns = l.length;
            double[] rhsValues = new double[this.numMetabs];
            char[] senses = new char[this.numMetabs];
            this.rxnExpressions = new GRBLinExpr[this.numMetabs];
            int met_count = false;

            int k;
            for(k = 0; k < this.numMetabs; ++k) {
                this.rxnExpressions[k] = new GRBLinExpr();
                senses[k] = '=';
                rhsValues[k] = 0.0;
            }

            for(k = 0; k < m.length; ++k) {
                Double cr = m[k][0];
                int cMet = cr.intValue();
                Double cc = m[k][1];
                int cVar = cc.intValue();
                this.rxnExpressions[cMet - 1].addTerm(m[k][2], this.rxnFluxes[cVar - 1]);
            }

            this.model.addConstrs(this.rxnExpressions, senses, rhsValues, (String[])null);
            this.model.update();
            this.setObjectiveReaction(this.numRxns, objIdxs, objWeights);
        } catch (GRBException var17) {
            System.out.println("Error in FBAOptimizerGurobi public constructor method");
            System.out.println("Error code: " + var17.getErrorCode() + ". " + var17.getMessage());
        }

    }

    public FBAOptimizerGurobi(double[][] m, double[] l, double[] u, int r) {
        this(m, l, u, new int[]{r});
    }

    public FBAOptimizerGurobi(double[][] m, double[] l, double[] u, int[] objIdxs) {
        this(m, l, u, objIdxs, getWeightsList(objIdxs));
    }

    protected static double[] getWeightsList(int[] idxs) {
        double[] weights = new double[idxs.length];

        for(int i = 0; i < weights.length; ++i) {
            weights[i] = idxs[i] >= 0 ? 1.0 : -1.0;
        }

        return weights;
    }

    public void initializeAbsModel(double[][] m, double[] l, double[] u, int[] objIdxs, double[] objWeights) {
        this.nVars = l.length * 2;
        this.createEmptyModelMin();
        this.createVarsModelMin(l, u);
        this.createObjFuncModelMin();
        this.addOrigConstraintsModelMin(m);
        this.addObjectiveListConstraintModelMin(objIdxs);
        this.addAbsSumConstraintsModelMin();
    }

    private void createEmptyModelMin() {
        try {
            this.envMin = new GRBEnv();
            this.envMin.set(IntParam.OutputFlag, 0);
            this.modelMin = new GRBModel(this.envMin);
        } catch (GRBException var2) {
            System.out.println("Error in FBAOptimizerGurobi.createEmptyModelMin");
            System.out.println("Error code: " + var2.getErrorCode() + ". " + var2.getMessage());
        }

    }

    private void createVarsModelMin(double[] l, double[] u) {
        double[] lbMin = new double[this.nVars];
        double[] ubMin = new double[this.nVars];
        double[] objMin = new double[this.nVars];
        char[] cTypeMin = new char[this.nVars];

        int k;
        for(k = 0; k < this.nVars / 2; ++k) {
            lbMin[k] = l[k];
            ubMin[k] = u[k];
            objMin[k] = 0.0;
            cTypeMin[k] = 'C';
        }

        for(k = this.nVars / 2; k < this.nVars; ++k) {
            lbMin[k] = 0.0;
            ubMin[k] = Double.MAX_VALUE;
            objMin[k] = 1.0;
            cTypeMin[k] = 'C';
        }

        try {
            this.modelMinVars = this.modelMin.addVars(lbMin, ubMin, objMin, cTypeMin, (String[])null);
            this.modelMin.update();
        } catch (GRBException var8) {
            System.out.println("Error in FBAOptimizerGurobi.createVarsModelMin");
            System.out.println("Error code: " + var8.getErrorCode() + ". " + var8.getMessage());
        }

    }

    private void createObjFuncModelMin() {
        GRBLinExpr objectiveFunc = new GRBLinExpr();

        for(int k = this.nVars / 2; k < this.nVars; ++k) {
            objectiveFunc.addTerm(1.0, this.modelMinVars[k]);
        }

        try {
            this.modelMin.setObjective(objectiveFunc, 1);
        } catch (GRBException var3) {
            System.out.println("Error in FBAOptimizerGurobi.createObjFuncModelMin");
            System.out.println("Error code: " + var3.getErrorCode() + ". " + var3.getMessage());
        }

    }

    private void sortByColumn(double[][] arr, int col) {
        Arrays.sort(arr, new Comparator<double[]>() {
            public int compare(double[] a, double[] b) {
                return Double.compare(a[0], b[0]);
            }
        });
    }

    private void addOrigConstraintsModelMin(double[][] m) {
        Double mtb = m[m.length - 1][0];
        int nMetabolites = mtb.intValue();
        GRBLinExpr[] origConstraints = new GRBLinExpr[nMetabolites];
        char[] senses = new char[nMetabolites];
        double[] rhs = new double[nMetabolites];

        int k;
        for(k = 0; k < nMetabolites; ++k) {
            origConstraints[k] = new GRBLinExpr();
            senses[k] = '=';
            rhs[k] = 0.0;
        }

        for(k = 0; k < m.length; ++k) {
            Double cr = m[k][0];
            int cMet = cr.intValue();
            Double cc = m[k][1];
            int cVar = cc.intValue();
            origConstraints[cMet - 1].addTerm(m[k][2], this.modelMinVars[cVar - 1]);
        }

        try {
            this.modelMin.addConstrs(origConstraints, senses, rhs, (String[])null);
            this.modelMin.update();
        } catch (GRBException var12) {
            System.out.println("Error in FBAOptimizerGurobi.addOrigConstraintsModelMin");
            System.out.println("Error code: " + var12.getErrorCode() + ". " + var12.getMessage());
        }

    }

    /** @deprecated */
    @Deprecated
    private void addBiomassConstraintModelMin(int biomassVarNum) {
        GRBLinExpr biomassConstraint = new GRBLinExpr();
        char biomassSense = '=';
        double biomassRHS = 0.0;
        double biomassCoef = 1.0;
        biomassConstraint.addTerm(biomassCoef, this.modelMinVars[biomassVarNum - 1]);

        try {
            this.modelMin.addConstr(biomassConstraint, biomassSense, biomassRHS, "biomassConstraint");
            this.modelMin.update();
        } catch (GRBException var9) {
            System.out.println("Error in FBAOptimizerGurobi.addBiomassConstraintModelMin");
            System.out.println("Error code: " + var9.getErrorCode() + ". " + var9.getMessage());
        }

    }

    private void addObjectiveListConstraintModelMin(int[] objIdxs) {
        int nObjs = objIdxs.length;
        char sense = '=';
        double rhs = 0.0;
        double coef = 1.0;

        for(int i = 0; i < nObjs; ++i) {
            String constraintName = "Obj" + String.valueOf(i);
            GRBLinExpr objConstraint = new GRBLinExpr();
            objConstraint.addTerm(coef, this.modelMinVars[objIdxs[i] - 1]);

            try {
                this.modelMin.addConstr(objConstraint, sense, rhs, constraintName);
                this.modelMin.update();
            } catch (GRBException var12) {
                System.out.println("Error in FBAOptimizerGurobi.addBiomassConstraintModelMin");
                System.out.println("Error code: " + var12.getErrorCode() + ". " + var12.getMessage());
            }
        }

    }

    private void addAbsSumConstraintsModelMin() {
        GRBLinExpr[] absSumConstraints = new GRBLinExpr[this.nVars];
        int counter = 0;
        int nOrigReactions = this.nVars / 2;
        char[] senses = new char[this.nVars];
        double[] rhs = new double[this.nVars];

        int k;
        for(k = 0; k < nOrigReactions; ++k) {
            absSumConstraints[counter] = new GRBLinExpr();
            absSumConstraints[counter].addTerm(1.0, this.modelMinVars[k]);
            absSumConstraints[counter].addTerm(-1.0, this.modelMinVars[k + nOrigReactions]);
            ++counter;
            absSumConstraints[counter] = new GRBLinExpr();
            absSumConstraints[counter].addTerm(-1.0, this.modelMinVars[k]);
            absSumConstraints[counter].addTerm(-1.0, this.modelMinVars[k + nOrigReactions]);
            ++counter;
        }

        for(k = 0; k < this.nVars; ++k) {
            senses[k] = '<';
            rhs[k] = 0.0;
        }

        try {
            this.modelMin.addConstrs(absSumConstraints, senses, rhs, (String[])null);
            this.modelMin.update();
        } catch (GRBException var7) {
            System.out.println("Error in FBAOptimizerGurobi.addAbsSumConstraintsModelMin");
            System.out.println("Error code: " + var7.getErrorCode() + ". " + var7.getMessage());
        }

    }

    public int setExchLowerBounds(int[] exch, double[] lb) {
        if (lb.length != exch.length) {
            return 0;
        } else {
            GRBVar[] exchFluxes = new GRBVar[exch.length];

            for(int i = 0; i < exch.length; ++i) {
                try {
                    this.rxnFluxes[exch[i] - 1].set(DoubleAttr.LB, lb[i]);
                    exchFluxes[i] = this.rxnFluxes[exch[i] - 1];
                } catch (GRBException var7) {
                    System.out.println("Error in FBAOptimizerGurobi.setExchLowerBounds, exch loop");
                    System.out.println("Error code: " + var7.getErrorCode() + ". " + var7.getMessage());
                }
            }

            try {
                this.model.set(DoubleAttr.LB, exchFluxes, lb);
                this.model.update();
            } catch (GRBException var6) {
                System.out.println("Error in FBAOptimizerGurobi.setExchLowerBounds, update step");
                System.out.println("Error code: " + var6.getErrorCode() + ". " + var6.getMessage());
            }

            this.setExchLowerBoundsModelMin(exch, lb);
            return 1;
        }
    }

    private void setExchLowerBoundsModelMin(int[] exch, double[] lb) {
        GRBVar[] exchFluxes = new GRBVar[exch.length];

        for(int k = 0; k < exch.length; ++k) {
            exchFluxes[k] = this.modelMinVars[exch[k] - 1];
        }

        try {
            this.modelMin.set(DoubleAttr.LB, exchFluxes, lb);
        } catch (GRBException var5) {
            System.out.println("Error in FBAOptimizerGurobi.setExchLowerBoundsModelMin");
            System.out.println("Error code: " + var5.getErrorCode() + ". " + var5.getMessage());
        }

    }

    public int setExchUpperBounds(int[] exch, double[] ub) {
        if (ub.length != exch.length) {
            return 0;
        } else {
            GRBVar[] exchFluxes = new GRBVar[exch.length];

            for(int i = 0; i < exch.length; ++i) {
                try {
                    this.rxnFluxes[exch[i]].set(DoubleAttr.UB, ub[i]);
                    exchFluxes[i] = this.rxnFluxes[exch[i] - 1];
                } catch (GRBException var7) {
                    System.out.println("Error in FBAOptimizerGurobi.setExchUpperBounds, exch loop");
                    System.out.println("Error code: " + var7.getErrorCode() + ". " + var7.getMessage());
                }
            }

            try {
                this.model.set(DoubleAttr.UB, exchFluxes, ub);
                this.model.update();
            } catch (GRBException var6) {
                System.out.println("Error in FBAOptimizerGurobi.setExchUpperBounds, update step");
                System.out.println("Error code: " + var6.getErrorCode() + ". " + var6.getMessage());
            }

            return 1;
        }
    }

    public int setLowerBounds(int nrxns, double[] lb) {
        if (nrxns == 0) {
            return 2;
        } else {
            for(int i = 0; i < nrxns; ++i) {
                try {
                    this.rxnFluxes[i].set(DoubleAttr.LB, lb[i]);
                } catch (GRBException var5) {
                    System.out.println("Error in FBAOptimizerGurobi.setLowerBounds");
                    System.out.println("Error code: " + var5.getErrorCode() + ". " + var5.getMessage());
                }
            }

            this.setLowerBoundsModelMin(nrxns, lb);
            return 1;
        }
    }

    public int setLowerBoundsModelMin(int nrxns, double[] lb) {
        if (nrxns == 0) {
            return 2;
        } else {
            for(int i = 0; i < nrxns; ++i) {
                try {
                    this.modelMinVars[i].set(DoubleAttr.LB, lb[i]);
                } catch (GRBException var5) {
                    System.out.println("Error in FBAOptimizerGurobi.setLowerBoundsModelMin");
                    System.out.println("Error code: " + var5.getErrorCode() + ". " + var5.getMessage());
                }
            }

            return 1;
        }
    }

    public double[] getLowerBounds(int nrxns) {
        GRBVar[] rxnFluxesLocal = this.model.getVars();
        double[] l = new double[nrxns];

        for(int i = 0; i < nrxns; ++i) {
            try {
                l[i] = rxnFluxesLocal[i].get(DoubleAttr.LB);
            } catch (GRBException var6) {
                System.out.println("Error in FBAOptimizerGurobi.getLowerBounds");
                System.out.println("Error code: " + var6.getErrorCode() + ". " + var6.getMessage());
            }
        }

        return l;
    }

    public int setUpperBounds(int nrxns, double[] ub) {
        if (nrxns == 0) {
            return 2;
        } else {
            for(int i = 0; i < nrxns; ++i) {
                try {
                    this.rxnFluxes[i].set(DoubleAttr.UB, ub[i]);
                } catch (GRBException var5) {
                    System.out.println("Error in FBAOptimizerGurobi.setUpperBounds");
                    System.out.println("Error code: " + var5.getErrorCode() + ". " + var5.getMessage());
                }
            }

            this.setUpperBoundsModelMin(nrxns, ub);
            return 1;
        }
    }

    public int setUpperBoundsModelMin(int nrxns, double[] ub) {
        if (nrxns == 0) {
            return 2;
        } else {
            for(int i = 0; i < nrxns; ++i) {
                try {
                    this.modelMinVars[i].set(DoubleAttr.UB, ub[i]);
                } catch (GRBException var5) {
                    System.out.println("Error in FBAOptimizerGurobi.setUpperBoundsModelMin");
                    System.out.println("Error code: " + var5.getErrorCode() + ". " + var5.getMessage());
                }
            }

            return 1;
        }
    }

    public double[] getUpperBounds(int nrxns) {
        GRBVar[] rxnFluxesLocal = this.model.getVars();
        double[] u = new double[nrxns];

        for(int i = 0; i < nrxns; ++i) {
            try {
                u[i] = rxnFluxesLocal[i].get(DoubleAttr.UB);
            } catch (GRBException var6) {
                System.out.println("Error in FBAOptimizerGurobi.getUpperBounds");
                System.out.println("Error code: " + var6.getErrorCode() + ". " + var6.getMessage());
            }
        }

        return u;
    }

    public int setObjectiveReaction(int nrxns, int r) {
        return this.setObjectiveReaction(nrxns, new int[]{r}, new double[]{1.0});
    }

    public int setObjectiveReaction(int numRxns, int[] objs) {
        double[] max = new double[objs.length];

        for(int i = 0; i < objs.length; ++i) {
            max[i] = objs[i] >= 0 ? 1.0 : -1.0;
        }

        return this.setObjectiveReaction(numRxns, objs, max);
    }

    public int setObjectiveReaction(int nrxns, int[] idxs, double[] weights) {
        int[] var4 = idxs;
        int var5 = idxs.length;

        int i;
        for(int var6 = 0; var6 < var5; ++var6) {
            i = var4[var6];
            if (i < 1 || i > nrxns) {
                return 0;
            }
        }

        if (this.numMetabs != 0 && this.numRxns != 0) {
            int nObjs = idxs.length;
            String[] objRxnNames = new String[nObjs];
            int[] priorities = new int[nObjs];

            for(i = 0; i < nObjs; ++i) {
                priorities[i] = nObjs - i;
            }

            try {
                this.model.set(IntAttr.NumObj, nObjs);
                this.model.set(IntAttr.ModelSense, -1);

                for(i = 0; i < nObjs; ++i) {
                    this.model.set(IntParam.ObjNumber, i);
                    this.model.set(DoubleAttr.ObjNWeight, weights[i]);
                    String vname = "Obj" + String.valueOf(i);
                    this.model.set(StringAttr.ObjNName, vname);
                    objRxnNames[i] = vname;
                    GRBLinExpr expr = new GRBLinExpr();
                    expr.addTerm(1.0, this.rxnFluxes[idxs[i] - 1]);
                    this.model.setObjectiveN(expr, i, 0, weights[i], 0.0, 0.0, vname);
                }
            } catch (GRBException var10) {
                System.out.println("Error in FBAOptimizerGurobi.setObjectiveReaction");
                System.out.println("Error code: " + var10.getErrorCode() + ". " + var10.getMessage());
            }

            this.objReactions = idxs;
            this.objConstraintNames = objRxnNames;
            this.objWeights = weights;
            return 1;
        } else {
            return 2;
        }
    }

    public synchronized int run(int objSty) {
        this.runSuccess = false;
        if (this.numMetabs != 0 && this.numRxns != 0) {
            byte ret;
            int status;
            ret = -1;
            status = -1;
            int optimstatus;
            int optimstatus_minflux;
            label72:
            switch (objSty) {
                case 0:
                    try {
                        this.setObjectiveReaction(this.numRxns, this.objReactions, this.objWeights);
                        this.model.update();
                        this.model.optimize();
                        status = this.model.get(IntAttr.Status);
                        this.rxnFluxes = this.model.getVars();
                        optimstatus = this.model.get(IntAttr.Status);
                        if (optimstatus == 2) {
                            for(optimstatus_minflux = 0; optimstatus_minflux < this.rxnFluxes.length; ++optimstatus_minflux) {
                                this.fluxesModel[optimstatus_minflux] = this.rxnFluxes[optimstatus_minflux].get(DoubleAttr.X);
                            }

                            ret = 0;
                        } else {
                            optimstatus_minflux = 0;

                            while(true) {
                                if (optimstatus_minflux >= this.rxnFluxes.length) {
                                    break label72;
                                }

                                this.fluxesModel[optimstatus_minflux] = 0.0;
                                ++optimstatus_minflux;
                            }
                        }
                    } catch (GRBException var9) {
                        System.out.println("Error in FBAOptimizerGurobi.run, case MAXIMIZE_OBJECTIVE_FLUX");
                        System.out.println("Error code: " + var9.getErrorCode() + ". " + var9.getMessage());
                    }
                    break;
                case 4:
                    try {
                        this.setObjectiveReaction(this.numRxns, this.objReactions, this.objWeights);
                        this.model.update();
                        this.model.optimize();
                        this.rxnFluxes = this.model.getVars();
                        optimstatus = this.model.get(IntAttr.Status);
                        if (optimstatus == 2) {
                            for(optimstatus_minflux = 0; optimstatus_minflux < this.objReactions.length; ++optimstatus_minflux) {
                                double maximizedObjective = this.rxnFluxes[this.objReactions[optimstatus_minflux] - 1].get(DoubleAttr.X);
                                this.setObjectiveFluxToSpecificValue(optimstatus_minflux, maximizedObjective);
                            }

                            this.modelMin.optimize();
                            status = this.model.get(IntAttr.Status);
                            optimstatus_minflux = this.model.get(IntAttr.Status);
                            if (optimstatus_minflux == 2) {
                                this.modelMinVars = this.modelMin.getVars();

                                for(int k = 0; k < this.modelMinVars.length / 2; ++k) {
                                    this.fluxesModel[k] = this.modelMinVars[k].get(DoubleAttr.X);
                                }

                                ret = 0;
                            }
                        }
                    } catch (GRBException var8) {
                        System.out.println("Error in FBAOptimizerGurobi.run, case MAX_OBJECTIVE_MIN_TOTAL ");
                        System.out.println("Error code: " + var8.getErrorCode() + ". " + var8.getMessage());
                    }
            }

            if (ret == 0) {
                this.runSuccess = true;
            }

            if (status == 2) {
                status = 5;
            } else if (status == 3) {
                status = 4;
            }

            return status;
        } else {
            return 2;
        }
    }

    private void setObjectiveFluxToSpecificValue(int objIdx, double objectiveFlux) {
        GRBConstr[] objFluxVar = new GRBConstr[1];
        int var10000 = this.objReactions[objIdx];

        try {
            objFluxVar[0] = this.modelMin.getConstrByName(this.objConstraintNames[objIdx]);
        } catch (GRBException var9) {
            System.out.println("Error in FBAOptimizerGurobi.setObjectiveFluxToSpecificValue, attempting modelMin.getConstrByName");
            System.out.println("Error code: " + var9.getErrorCode() + ". " + var9.getMessage());
        }

        double[] v = new double[]{objectiveFlux};

        try {
            this.modelMin.set(DoubleAttr.RHS, objFluxVar, v);
        } catch (GRBException var8) {
            System.out.println("   Error in setObjectiveFluxToSpecificValue");
            System.out.println("Error in FBAOptimizerGurobi.setObjectiveFluxToSpecificValue, attempting to modelMin.set");
            System.out.println("Error code: " + var8.getErrorCode() + ". " + var8.getMessage());
        }

    }

    public double[] getFluxes() {
        double[] v = new double[this.numRxns];
        if (this.runSuccess) {
            for(int i = 0; i < this.numRxns; ++i) {
                v[i] = this.fluxesModel[i];
            }
        }

        return v;
    }

    public double[] getExchangeFluxes(int[] exch) {
        double[] v = new double[exch.length];
        if (this.runSuccess) {
            for(int i = 0; i < exch.length; ++i) {
                v[i] = this.fluxesModel[exch[i] - 1];
            }
        }

        return v;
    }

    public double getObjectiveSolution(int objreact) {
        try {
            if (this.runSuccess) {
                return this.rxnFluxes[objreact - 1].get(DoubleAttr.X);
            }
        } catch (GRBException var3) {
            System.out.println("Error in FBAOptimizerGurobi.getObjectiveSolution");
            System.out.println("Error code: " + var3.getErrorCode() + ". " + var3.getMessage());
        }

        return -1.7976931348623157E308;
    }

    public double[] getObjectiveSolutions(int[] objs) {
        double[] res = new double[objs.length];

        for(int i = 0; i < objs.length; ++i) {
            res[i] = this.getObjectiveSolution(objs[i]);
        }

        return res;
    }

    public double getObjectiveFluxSolution(int objreact) {
        try {
            if (this.runSuccess) {
                return this.rxnFluxes[objreact - 1].get(DoubleAttr.X);
            }
        } catch (GRBException var3) {
            System.out.println("Error in FBAOptimizerGurobi.getObjectiveFluxSolution");
            System.out.println("Error code: " + var3.getErrorCode() + ". " + var3.getMessage());
        }

        return -1.7976931348623157E308;
    }

    public int getFBAstatus() {
        return this.runSuccess ? 1 : 0;
    }

    public int setObjectiveUpperBound(int objreact, double ub) {
        GRBVar[] objFlux = new GRBVar[1];
        double[] ubarray = new double[1];

        try {
            this.rxnFluxes[objreact - 1].set(DoubleAttr.UB, ub);
            objFlux[0] = this.rxnFluxes[objreact - 1];
            ubarray[0] = ub;
            this.model.set(DoubleAttr.UB, objFlux, ubarray);
            this.model.update();
        } catch (GRBException var7) {
            System.out.println("Error in FBAOptimizerGurobi.setObjectiveUpperBound");
            System.out.println("Error code: " + var7.getErrorCode() + ". " + var7.getMessage());
        }

        return 1;
    }

    public int setObjectiveLowerBound(int objreact, double lb) {
        GRBVar[] objFlux = new GRBVar[1];
        double[] lbarray = new double[1];

        try {
            this.rxnFluxes[objreact - 1].set(DoubleAttr.LB, lb);
            objFlux[0] = this.rxnFluxes[objreact - 1];
            lbarray[0] = lb;
            this.model.set(DoubleAttr.LB, objFlux, lbarray);
            this.model.update();
        } catch (GRBException var7) {
            System.out.println("Error in FBAOptimizerGurobi.setObjectiveLowerBound");
            System.out.println("Error code: " + var7.getErrorCode() + ". " + var7.getMessage());
        }

        return 1;
    }

    public FBAOptimizerGurobi clone() {
        GRBVar[] rxnFluxesLocal = this.model.getVars();
        double[] l = new double[rxnFluxesLocal.length];
        double[] u = new double[rxnFluxesLocal.length];
        int[] objs = this.objReactions;
        double[] weights = this.objWeights;

        try {
            for(int k = 0; k < rxnFluxesLocal.length; ++k) {
                l[k] = rxnFluxesLocal[k].get(DoubleAttr.LB);
                u[k] = rxnFluxesLocal[k].get(DoubleAttr.UB);
            }
        } catch (GRBException var7) {
            System.out.println("Error in FBAOptimizerGurobi.clone");
            System.out.println("Error code: " + var7.getErrorCode() + ". " + var7.getMessage());
        }

        FBAOptimizerGurobi optimizerCopy = new FBAOptimizerGurobi(this.stoichMatrix, l, u, objs, weights);
        return optimizerCopy;
    }

    public int setObjectiveWeights(double[] objWt) {
        this.objWeights = objWt;
        return 1;
    }
}
