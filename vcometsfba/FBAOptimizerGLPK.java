//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package edu.bu.segrelab.comets.fba;

import edu.bu.segrelab.comets.CometsConstants;
import org.gnu.glpk.GLPK;
import org.gnu.glpk.GLPKConstants;
import org.gnu.glpk.SWIGTYPE_p_double;
import org.gnu.glpk.SWIGTYPE_p_f_p_void_p_q_const__char__int;
import org.gnu.glpk.SWIGTYPE_p_int;
import org.gnu.glpk.SWIGTYPE_p_void;
import org.gnu.glpk.glp_iptcp;
import org.gnu.glpk.glp_prob;
import org.gnu.glpk.glp_smcp;

public class FBAOptimizerGLPK extends FBAOptimizer implements CometsConstants {
    public static final int SIMPLEX_METHOD = 0;
    public static final int INTERIOR_POINT_METHOD = 1;
    private int glpkSolverMethod;
    private glp_prob lp;
    private glp_prob lpMSA;
    private glp_smcp simParam;
    private glp_iptcp intParam;
    private int numRxns;
    private int numMetabs;
    private int numExch;
    private boolean runSuccess;
    SWIGTYPE_p_int forceIdx;
    SWIGTYPE_p_double forceVal;
    SWIGTYPE_p_int lpRow;
    SWIGTYPE_p_int msaRow;
    private int objReaction;
    private int objStyle;
    private static int[] GLPIntParam;
    private static double[] GLPRealParam;

    private FBAOptimizerGLPK() {
        this.glpkSolverMethod = 0;
        this.runSuccess = false;
        this.lp = GLPK.glp_create_prob();
        this.lpMSA = GLPK.glp_create_prob();
        this.forceIdx = GLPK.new_intArray(2);
        this.forceVal = GLPK.new_doubleArray(2);
        this.lpRow = GLPK.new_intArray(1);
        this.msaRow = GLPK.new_intArray(1);
        GLPK.glp_term_hook((SWIGTYPE_p_f_p_void_p_q_const__char__int)null, (SWIGTYPE_p_void)null);
        GLPK.glp_set_obj_dir(this.lp, GLPK.GLP_MAX);
        GLPK.glp_set_obj_dir(this.lpMSA, GLPK.GLP_MAX);
        this.objStyle = 0;
        this.setParameters();
    }

    public FBAOptimizerGLPK(double[][] m, double[] l, double[] u, int r) {
        this();
        this.setStoichiometricMatrix(m);
        this.setBounds(l, u);
        this.setObjectiveReaction(this.numRxns, r);
    }

    public int getSolverMethod() {
        return this.glpkSolverMethod;
    }

    public void setSolverMethod(int method) {
        if (method == 0 || method == 1) {
            this.glpkSolverMethod = method;
        }
    }

    public void setModelName(String name) {
        GLPK.glp_set_prob_name(this.lp, name);
        GLPK.glp_set_prob_name(this.lpMSA, name + "_MSA");
    }

    public String getModelName() {
        return GLPK.glp_get_prob_name(this.lp);
    }

    private void setStoichiometricMatrix(double[][] m) {
        this.numMetabs = m.length;
        this.numRxns = m[0].length;
        GLPK.glp_add_cols(this.lp, this.numRxns);
        GLPK.glp_add_cols(this.lpMSA, this.numRxns * 2);

        int i;
        for(i = 0; i < this.numRxns; ++i) {
            GLPK.glp_set_col_name(this.lp, i + 1, "r" + (i + 1));
            GLPK.glp_set_col_kind(this.lp, i + 1, GLPK.GLP_CV);
            GLPK.glp_set_col_name(this.lpMSA, i + 1, "r" + (i + 1));
            GLPK.glp_set_col_kind(this.lpMSA, i + 1, GLPK.GLP_CV);
            GLPK.glp_set_col_name(this.lpMSA, i + this.numRxns + 1, "d" + (i + 1));
            GLPK.glp_set_col_kind(this.lpMSA, i + this.numRxns + 1, GLPK.GLP_CV);
        }

        GLPK.glp_add_rows(this.lp, this.numMetabs);
        GLPK.glp_add_rows(this.lpMSA, this.numMetabs + this.numRxns + this.numRxns);

        for(i = 0; i < this.numMetabs; ++i) {
            GLPK.glp_set_row_bnds(this.lp, i + 1, GLPK.GLP_FX, 0.0, 0.0);
            GLPK.glp_set_row_bnds(this.lpMSA, i + 1, GLPK.GLP_FX, 0.0, 0.0);
        }

        for(i = this.numMetabs; i < this.numMetabs + 2 * this.numRxns; ++i) {
            GLPK.glp_set_row_bnds(this.lpMSA, i + 1, GLPK.GLP_UP, 0.0, 0.0);
        }

        SWIGTYPE_p_int ia = GLPK.new_intArray(this.numMetabs * this.numRxns + 1);
        SWIGTYPE_p_int ja = GLPK.new_intArray(this.numMetabs * this.numRxns + 1);
        SWIGTYPE_p_double mLin = GLPK.new_doubleArray(this.numMetabs * this.numRxns + 1);
        int ne = 0;

        for(int i = 0; i < this.numMetabs; ++i) {
            for(int j = 0; j < this.numRxns; ++j) {
                if (m[i][j] != 0.0) {
                    ++ne;
                    GLPK.intArray_setitem(ia, ne, i + 1);
                    GLPK.intArray_setitem(ja, ne, j + 1);
                    GLPK.doubleArray_setitem(mLin, ne, m[i][j]);
                }
            }
        }

        GLPK.glp_load_matrix(this.lp, ne, ia, ja, mLin);
        SWIGTYPE_p_int iaMSA = GLPK.new_intArray(this.numMetabs * this.numRxns + 4 * this.numRxns + 1);
        SWIGTYPE_p_int jaMSA = GLPK.new_intArray(this.numMetabs * this.numRxns + 4 * this.numRxns + 1);
        SWIGTYPE_p_double mLinMSA = GLPK.new_doubleArray(this.numMetabs * this.numRxns + 4 * this.numRxns + 1);
        ne = 0;

        int i;
        for(i = 0; i < this.numMetabs + 2 * this.numRxns; ++i) {
            for(int j = 0; j < 2 * this.numRxns; ++j) {
                if (i < this.numMetabs && j < this.numRxns && m[i][j] != 0.0) {
                    ++ne;
                    GLPK.intArray_setitem(iaMSA, ne, i + 1);
                    GLPK.intArray_setitem(jaMSA, ne, j + 1);
                    GLPK.doubleArray_setitem(mLinMSA, ne, m[i][j]);
                } else if (i >= this.numMetabs || j < this.numRxns) {
                    if (i >= this.numMetabs && i < this.numMetabs + this.numRxns && j < this.numRxns && i - this.numMetabs == j) {
                        ++ne;
                        GLPK.intArray_setitem(iaMSA, ne, i + 1);
                        GLPK.intArray_setitem(jaMSA, ne, j + 1);
                        GLPK.doubleArray_setitem(mLinMSA, ne, 1.0);
                    } else if (i >= this.numMetabs && i < this.numMetabs + this.numRxns && j >= this.numRxns && i - this.numMetabs == j - this.numRxns) {
                        ++ne;
                        GLPK.intArray_setitem(iaMSA, ne, i + 1);
                        GLPK.intArray_setitem(jaMSA, ne, j + 1);
                        GLPK.doubleArray_setitem(mLinMSA, ne, -1.0);
                    } else if (i >= this.numMetabs + this.numRxns && (i - this.numMetabs - this.numRxns == j || i - this.numMetabs - this.numRxns == j - this.numRxns)) {
                        ++ne;
                        GLPK.intArray_setitem(iaMSA, ne, i + 1);
                        GLPK.intArray_setitem(jaMSA, ne, j + 1);
                        GLPK.doubleArray_setitem(mLinMSA, ne, -1.0);
                    }
                }
            }
        }

        GLPK.glp_load_matrix(this.lpMSA, ne, iaMSA, jaMSA, mLinMSA);

        for(i = 1; i <= this.numRxns * 2; ++i) {
            GLPK.glp_set_obj_coef(this.lpMSA, i, 0.0);
        }

        for(i = this.numRxns + 1; i <= this.numRxns * 2; ++i) {
            GLPK.glp_set_col_bnds(this.lpMSA, i, GLPK.GLP_DB, 0.0, Double.MAX_VALUE);
            GLPK.glp_set_obj_coef(this.lpMSA, i, 1.0);
        }

        GLPK.glp_set_obj_name(this.lpMSA, "sum_values");
    }

    public int setBounds(double[] lb, double[] ub) {
        if (lb.length != ub.length) {
            return 0;
        } else if (this.numMetabs != 0 && this.numRxns != 0) {
            int i;
            for(i = 0; i < lb.length; ++i) {
                if (lb[i] > ub[i]) {
                    return 0;
                }
            }

            for(i = 0; i < lb.length; ++i) {
                int type = GLPKConstants.GLP_DB;
                if (lb[i] == ub[i]) {
                    type = GLPKConstants.GLP_FX;
                }

                GLPK.glp_set_col_bnds(this.lp, i + 1, type, lb[i], ub[i]);
                GLPK.glp_set_col_bnds(this.lpMSA, i + 1, type, lb[i], ub[i]);
            }

            return 1;
        } else {
            return 2;
        }
    }

    public int setLowerBounds(int nrxns, double[] lb) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            if (lb.length != nrxns) {
                return 0;
            } else {
                for(int i = 0; i < this.numRxns; ++i) {
                    double u = GLPK.glp_get_col_ub(this.lp, i + 1);
                    int type = GLPKConstants.GLP_DB;
                    if (lb[i] == u) {
                        type = GLPKConstants.GLP_FX;
                    }

                    GLPK.glp_set_col_bnds(this.lp, i + 1, type, lb[i], u);
                    GLPK.glp_set_col_bnds(this.lpMSA, i + 1, type, lb[i], u);
                }

                return 1;
            }
        } else {
            return 2;
        }
    }

    public int setExchLowerBounds(int[] exch, double[] lb) {
        for(int i = 0; i < exch.length; ++i) {
            double u = GLPK.glp_get_col_ub(this.lp, exch[i]);
            int type = GLPKConstants.GLP_DB;
            if (lb[i] == u) {
                type = GLPKConstants.GLP_FX;
            }

            GLPK.glp_set_col_bnds(this.lp, exch[i], type, lb[i], u);
            GLPK.glp_set_col_bnds(this.lpMSA, exch[i], type, lb[i], u);
        }

        return 1;
    }

    public int setExchUpperBounds(int[] exch, double[] ub) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            if (ub.length != this.numExch) {
                return 0;
            } else {
                for(int i = 0; i < this.numExch; ++i) {
                    double l = GLPK.glp_get_col_ub(this.lp, exch[i]);
                    int type = GLPKConstants.GLP_DB;
                    if (ub[i] == l) {
                        type = GLPKConstants.GLP_FX;
                    }

                    GLPK.glp_set_col_bnds(this.lp, exch[i], type, l, ub[i]);
                    GLPK.glp_set_col_bnds(this.lpMSA, exch[i], type, l, ub[i]);
                }

                return 1;
            }
        } else {
            return 2;
        }
    }

    public double[] getLowerBounds(int nrxns) {
        double[] l = new double[nrxns];

        for(int i = 0; i < nrxns; ++i) {
            l[i] = GLPK.glp_get_col_lb(this.lp, i + 1);
        }

        return l;
    }

    public int setUpperBounds(int nrxns, double[] ub) {
        if (this.numMetabs != 0 && this.numRxns != 0) {
            for(int i = 0; i < nrxns; ++i) {
                double l = GLPK.glp_get_col_lb(this.lp, i + 1);
                int type = GLPKConstants.GLP_DB;
                if (l == ub[i]) {
                    type = GLPKConstants.GLP_FX;
                }

                GLPK.glp_set_col_bnds(this.lp, i + 1, type, l, ub[i]);
                GLPK.glp_set_col_bnds(this.lpMSA, i + 1, type, l, ub[i]);
            }

            return 1;
        } else {
            return 2;
        }
    }

    public double[] getUpperBounds(int nrxns) {
        double[] u = new double[nrxns];

        for(int i = 0; i < nrxns; ++i) {
            u[i] = GLPK.glp_get_col_ub(this.lp, i + 1);
        }

        return u;
    }

    public int setObjectiveStyle(int obj) {
        if (obj != 0 && obj != 1 && obj != 2 && obj != 3 && obj != 4 && obj != 5 && obj != 6 && obj != 7) {
            return 0;
        } else {
            this.objStyle = obj;
            if (this.objStyle != 2 && this.objStyle != 3) {
                this.setObjectiveReaction(this.numRxns, this.objReaction);
            }

            return 1;
        }
    }

    public int getObjectiveStyle() {
        return this.objStyle;
    }

    public int getObjectiveIndex() {
        return this.objReaction;
    }

    public int setObjectiveReaction(int nrxns, int r) {
        if (r >= 1 && r <= nrxns) {
            if (this.numMetabs != 0 && nrxns != 0) {
                for(int i = 1; i <= nrxns; ++i) {
                    GLPK.glp_set_obj_coef(this.lp, i, 0.0);
                }

                GLPK.glp_set_obj_coef(this.lp, r, 1.0);
                GLPK.glp_set_obj_name(this.lp, "obj");
                this.objReaction = r;
                return 1;
            } else {
                return 2;
            }
        } else {
            return 0;
        }
    }

    public void delete() {
        GLPK.glp_delete_prob(this.lp);
        GLPK.glp_delete_prob(this.lpMSA);
    }

    public void setParameters() {
        this.simParam = new glp_smcp();
        GLPK.glp_init_smcp(this.simParam);
        this.intParam = new glp_iptcp();
        GLPK.glp_init_iptcp(this.intParam);
        GLPK.glp_init_smcp(this.simParam);
        GLPK.glp_init_iptcp(this.intParam);
        this.simParam.setMsg_lev(GLPKConstants.GLP_MSG_OFF);
        this.intParam.setMsg_lev(GLPKConstants.GLP_MSG_OFF);
        switch (GLPIntParam[2]) {
            case 0:
                this.simParam.setMeth(GLPKConstants.GLP_PRIMAL);
                break;
            case 1:
                this.simParam.setMeth(GLPKConstants.GLP_DUAL);
                break;
            case 2:
                this.simParam.setMeth(GLPKConstants.GLP_DUALP);
        }

        if (GLPIntParam[3] == 0) {
            this.simParam.setPricing(GLPKConstants.GLP_PT_STD);
        } else {
            this.simParam.setPricing(GLPKConstants.GLP_PT_PSE);
        }

        if (GLPIntParam[20] == 0) {
            this.simParam.setR_test(GLPKConstants.GLP_RT_STD);
        } else {
            this.simParam.setR_test(GLPKConstants.GLP_RT_HAR);
        }

        this.simParam.setTol_bnd(GLPRealParam[1]);
        this.simParam.setTol_dj(GLPRealParam[2]);
        this.simParam.setTol_piv(GLPRealParam[3]);
        this.simParam.setObj_ll(GLPRealParam[4]);
        this.simParam.setObj_ul(GLPRealParam[5]);
        if (GLPIntParam[5] == -1) {
            this.simParam.setIt_lim(Integer.MAX_VALUE);
        } else {
            this.simParam.setIt_lim(GLPIntParam[5]);
        }

        if (GLPRealParam[6] == -1.0) {
            this.simParam.setTm_lim(Integer.MAX_VALUE);
        } else {
            this.simParam.setTm_lim((int)GLPRealParam[6]);
        }

        this.simParam.setOut_frq(GLPIntParam[7]);
        this.simParam.setOut_dly((int)GLPRealParam[7]);
        if (GLPIntParam[16] != 0) {
            this.simParam.setPresolve(GLPK.GLP_ON);
        } else {
            this.simParam.setPresolve(GLPK.GLP_OFF);
        }

    }

    public synchronized int run(int objSty) {
        this.runSuccess = false;
        if (this.numMetabs != 0 && this.numRxns != 0) {
            int ret;
            this.objStyle = objSty;
            int ret = true;
            label34:
            switch (objSty) {
                case 0:
                    GLPK.glp_set_obj_dir(this.lp, GLPK.GLP_MAX);
                    switch (this.glpkSolverMethod) {
                        case 0:
                            ret = GLPK.glp_simplex(this.lp, this.simParam);
                            break label34;
                        default:
                            ret = GLPK.glp_interior(this.lp, this.intParam);
                            break label34;
                    }
                case 1:
                    GLPK.glp_set_obj_dir(this.lp, GLPK.GLP_MIN);
                    switch (this.glpkSolverMethod) {
                        case 0:
                            ret = GLPK.glp_simplex(this.lp, this.simParam);
                            break label34;
                        default:
                            ret = GLPK.glp_interior(this.lp, this.intParam);
                            break label34;
                    }
                case 2:
                    ret = this.runOptimizeSumFluxes();
                    break;
                case 3:
                    ret = this.runOptimizeSumFluxes();
                    break;
                default:
                    ret = this.runOptimizeSumAbsoluteValuesFluxes();
            }

            if (ret == 0) {
                this.runSuccess = true;
            }

            return objSty != 0 && objSty != 1 && objSty != 2 && objSty != 3 ? GLPK.glp_get_status(this.lpMSA) : GLPK.glp_get_status(this.lp);
        } else {
            return 2;
        }
    }

    private synchronized int runOptimizeSumFluxes() {
        this.setObjectiveReaction(this.numRxns, this.objReaction);
        GLPK.glp_set_obj_dir(this.lp, GLPK.GLP_MAX);
        int ret = true;
        int ret;
        switch (this.glpkSolverMethod) {
            case 0:
                ret = GLPK.glp_simplex(this.lp, this.simParam);
                break;
            default:
                ret = GLPK.glp_interior(this.lp, this.intParam);
        }

        if (ret != 0) {
            return ret;
        } else {
            double sol = GLPK.glp_get_obj_val(this.lp);
            double lb = GLPK.glp_get_col_lb(this.lp, this.objReaction);
            double ub = GLPK.glp_get_col_ub(this.lp, this.objReaction);
            GLPK.glp_set_col_bnds(this.lp, this.objReaction, GLPK.GLP_FX, sol, sol);

            int type;
            for(type = 1; type <= this.numRxns; ++type) {
                GLPK.glp_set_obj_coef(this.lp, type, 1.0);
            }

            if (this.objStyle == 2) {
                GLPK.glp_set_obj_dir(this.lp, GLPK.GLP_MAX);
            } else {
                GLPK.glp_set_obj_dir(this.lp, GLPK.GLP_MIN);
            }

            switch (this.glpkSolverMethod) {
                case 0:
                    ret = GLPK.glp_simplex(this.lp, this.simParam);
                    break;
                default:
                    ret = GLPK.glp_interior(this.lp, this.intParam);
            }

            type = GLPK.GLP_DB;
            if (lb == ub) {
                type = GLPK.GLP_FX;
            }

            GLPK.glp_set_col_bnds(this.lp, this.objReaction, type, lb, ub);
            return ret;
        }
    }

    private synchronized int runOptimizeSumAbsoluteValuesFluxes() {
        int ret = true;
        if (this.objStyle != 7 && this.objStyle != 6) {
            GLPK.glp_set_obj_dir(this.lp, GLPK.GLP_MAX);
        } else {
            GLPK.glp_set_obj_dir(this.lp, GLPK.GLP_MIN);
        }

        if (this.objStyle != 4 && this.objStyle != 6) {
            GLPK.glp_set_obj_dir(this.lpMSA, GLPK.GLP_MAX);
        } else {
            GLPK.glp_set_obj_dir(this.lpMSA, GLPK.GLP_MIN);
        }

        int status = true;
        this.simParam.setPresolve(GLPK.GLP_OFF);
        int status;
        switch (this.glpkSolverMethod) {
            case 0:
                status = GLPK.glp_simplex(this.lp, this.simParam);
                break;
            default:
                status = GLPK.glp_interior(this.lp, this.intParam);
        }

        if (status != 0) {
            return status;
        } else {
            double sol = GLPK.glp_get_obj_val(this.lp);
            GLPK.glp_set_col_bnds(this.lpMSA, this.objReaction, GLPK.GLP_FX, sol, sol);
            int ret;
            switch (this.glpkSolverMethod) {
                case 0:
                    ret = GLPK.glp_simplex(this.lpMSA, this.simParam);
                    break;
                default:
                    ret = GLPK.glp_interior(this.lpMSA, this.intParam);
            }

            return ret;
        }
    }

    public double[] getFluxes() {
        double[] v = new double[this.numRxns];
        if (this.runSuccess) {
            for(int i = 0; i < this.numRxns; ++i) {
                if (this.objStyle != 0 && this.objStyle != 1 && this.objStyle != 2 && this.objStyle != 3) {
                    v[i] = GLPK.glp_get_col_prim(this.lpMSA, i + 1);
                } else {
                    v[i] = GLPK.glp_get_col_prim(this.lp, i + 1);
                }
            }
        }

        return v;
    }

    public double[] getExchangeFluxes(int[] exch) {
        double[] v = new double[exch.length];
        if (this.runSuccess) {
            for(int i = 0; i < exch.length; ++i) {
                if (this.objStyle != 0 && this.objStyle != 1 && this.objStyle != 2 && this.objStyle != 3) {
                    v[i] = GLPK.glp_get_col_prim(this.lpMSA, exch[i]);
                } else {
                    v[i] = GLPK.glp_get_col_prim(this.lp, exch[i]);
                }
            }
        }

        return v;
    }

    public void setMediaConditions(double[] media, int[] exch) {
        if (media.length == this.numExch) {
            for(int i = 0; i < media.length; ++i) {
                GLPK.glp_set_col_bnds(this.lp, exch[i], GLPKConstants.GLP_DB, media[i], GLPK.glp_get_col_ub(this.lp, exch[i]));
            }

        }
    }

    public double getObjectiveSolution(int r) {
        if (this.runSuccess) {
            return this.objStyle != 0 && this.objStyle != 1 ? GLPK.glp_get_obj_val(this.lpMSA) : GLPK.glp_get_obj_val(this.lp);
        } else {
            return -1.7976931348623157E308;
        }
    }

    public double getObjectiveFluxSolution(int objreact) {
        return this.runSuccess ? GLPK.glp_get_col_prim(this.lp, objreact) : -1.7976931348623157E308;
    }

    public int setObjectiveUpperBound(int objreact, double ub) {
        double lb = GLPK.glp_get_col_lb(this.lp, this.objReaction);
        if (ub < lb) {
            return 0;
        } else {
            int type = GLPKConstants.GLP_DB;
            if (lb == ub) {
                type = GLPKConstants.GLP_FX;
            }

            GLPK.glp_set_col_bnds(this.lp, objreact, type, lb, ub);
            return 1;
        }
    }

    public int setObjectiveLowerBound(int objreact, double lb) {
        return 1;
    }

    protected FBAOptimizerGLPK(glp_prob lp, glp_prob lpMSA, int numMetabs, int numRxns, int numExch) {
        this.glpkSolverMethod = 0;
        this.lp = lp;
        this.lpMSA = lpMSA;
        this.numMetabs = numMetabs;
        this.numRxns = numRxns;
        this.numExch = numExch;
    }

    public FBAOptimizerGLPK(double[][] m, double[] l, double[] u, int[] objs) {
        this(m, l, u, objs[0]);
        if (objs.length > 1) {
            System.out.println("Warning: the GLPK optimizer is not configured to support multiple objective reactions!");
            System.out.println("Using reaction " + String.valueOf(objs[0]) + " as the objective.");
        }

    }

    public FBAOptimizerGLPK clone() {
        glp_prob lpNew = GLPK.glp_create_prob();
        GLPK.glp_copy_prob(lpNew, this.lp, GLPK.GLP_ON);
        glp_prob lpMSANew = GLPK.glp_create_prob();
        GLPK.glp_copy_prob(lpMSANew, this.lpMSA, GLPK.GLP_ON);
        FBAOptimizerGLPK optimizerCopy = new FBAOptimizerGLPK(lpNew, lpMSANew, this.numMetabs, this.numRxns, this.numExch);
        optimizerCopy.setParameters();
        return optimizerCopy;
    }

    public int getFBAstatus() {
        return this.runSuccess ? 1 : 0;
    }

    public int setObjectiveReaction(int numRxns, int[] objs) {
        if (objs.length > 1) {
            System.out.println("Warning: the GLPK optimizer is not configured to support multiple objective reactions!");
            System.out.println("Using reaction " + String.valueOf(objs[0]) + " as the objective.");
        }

        return this.setObjectiveReaction(numRxns, objs[0]);
    }

    public double[] getObjectiveSolutions(int[] objReactions) {
        double[] res = new double[objReactions.length];

        for(int i = 0; i < objReactions.length; ++i) {
            res[i] = this.getObjectiveSolution(objReactions[i]);
        }

        return res;
    }

    public int setObjectiveWeights(double[] objWt) {
        return 0;
    }

    static {
        try {
            System.loadLibrary("glpk_java");
        } catch (UnsatisfiedLinkError var3) {
            UnsatisfiedLinkError e = var3;

            try {
                System.out.println(e);
                System.loadLibrary("glpk_4_47");
                System.loadLibrary("glpk_4_47_java");
            } catch (UnsatisfiedLinkError var2) {
                System.out.println(var2);
                System.loadLibrary("glpk_4_44");
                System.loadLibrary("glpk_4_44_java");
            }
        }

        GLPIntParam = new int[]{0, 1, 0, 1, 0, Integer.MAX_VALUE, Integer.MAX_VALUE, 200, 1, 2, 0, 1, 0, 0, 3, 2, 1, 0, 2, 0, 1};
        GLPRealParam = new double[]{0.07, 1.0E-7, 1.0E-7, 1.0E-10, -1.7976931348623157E308, Double.MAX_VALUE, 2.147483647E9, 0.0, 1.0E-5, 1.0E-7, 0.0};
    }
}
