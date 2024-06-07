//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package edu.bu.segrelab.comets.fba;

import edu.bu.segrelab.comets.Cell;
import edu.bu.segrelab.comets.Comets;
import edu.bu.segrelab.comets.CometsConstants;
import edu.bu.segrelab.comets.CometsParameters;
import edu.bu.segrelab.comets.Model;
import edu.bu.segrelab.comets.PackageParameters;
import edu.bu.segrelab.comets.World2D;
import edu.bu.segrelab.comets.World3D;
import edu.bu.segrelab.comets.util.Utility;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import jdistlib.Gamma;
import jdistlib.Poisson;
import jdistlib.rng.MersenneTwister;
import org.apache.commons.lang3.ArrayUtils;

public class FBACell extends Cell implements CometsConstants {
    private FBAModel[] fbaModels;
    private FBAWorld world;
    private FBAWorld3D world3D;
    private final int id;
    private int cellColor;
    private double[] biomass;
    private double[] convectionRHS1;
    private double[] convectionRHS2;
    private double jointRHS1;
    private double jointRHS2;
    private double[] deltaBiomass;
    private double[] dyingBiomass;
    private double[] allModelsGrowthRates;
    private double[][] deltaMedia;
    private boolean stationaryStatus;
    private double[][] fluxes;
    private int[] FBAstatus;
    private CometsParameters cParams;
    private FBAParameters pParams;
    private PrintWriter PoissWriter;

    public FBACell(int x, int y, FBAWorld world, FBAModel[] fbaModels, CometsParameters cParams, FBAParameters pParams) {
        this(x, y, new double[fbaModels.length], world, fbaModels, cParams, pParams);

        for(int i = 0; i < this.biomass.length; ++i) {
            this.biomass[i] = Utility.randomDouble();
        }

    }

    public FBACell(int x, int y, double[] biomass, FBAWorld world, FBAModel[] fbaModels, CometsParameters cParams, FBAParameters pParams) {
        this.stationaryStatus = false;
        this.x = x;
        this.y = y;
        this.id = getNewCellID();
        this.deltaBiomass = new double[biomass.length];
        this.dyingBiomass = new double[biomass.length];
        this.deltaMedia = new double[biomass.length][];
        this.FBAstatus = new int[biomass.length];
        this.fbaModels = fbaModels;
        this.world = world;
        this.cParams = cParams;
        this.pParams = pParams;
        this.biomass = new double[fbaModels.length];
        this.setBiomass(biomass);
        this.convectionRHS1 = new double[biomass.length];
        this.convectionRHS2 = new double[biomass.length];
        if (cParams.showGraphics()) {
            this.cellColor = this.calculateColor();
        } else {
            this.cellColor = 0;
        }

        this.fluxes = new double[fbaModels.length][];
        world.putCellAt(x, y, this);
        this.updateDiffusibility();
    }

    public FBACell(int x, int y, int z, double[] biomass, FBAWorld3D world3D, FBAModel[] fbaModels, CometsParameters cParams, FBAParameters pParams) {
        this.stationaryStatus = false;
        this.x = x;
        this.y = y;
        this.z = z;
        this.id = getNewCellID();
        this.deltaBiomass = new double[biomass.length];
        this.dyingBiomass = new double[biomass.length];
        this.FBAstatus = new int[biomass.length];
        this.fbaModels = fbaModels;
        this.world3D = world3D;
        this.cParams = cParams;
        this.pParams = pParams;
        this.biomass = new double[fbaModels.length];
        this.setBiomass3D(biomass);
        this.convectionRHS1 = new double[biomass.length];
        this.convectionRHS2 = new double[biomass.length];
        if (cParams.showGraphics()) {
            this.cellColor = this.calculateColor();
        } else {
            this.cellColor = 0;
        }

        this.fluxes = new double[fbaModels.length][];
        world3D.putCellAt(x, y, z, this);
        this.updateDiffusibility3D();
    }

    public int changeBiomass(double[] delta) {
        return this.setBiomass(delta, true);
    }

    public int setBiomass(double[] values) {
        return this.setBiomass(values, false);
    }

    public int setBiomass3D(double[] values) {
        return this.setBiomass3D(values, false);
    }

    public void setConvectionRHS1(double[] values) {
        for(int i = 0; i < this.convectionRHS1.length; ++i) {
            this.convectionRHS1[i] = values[i];
        }

    }

    public void setConvectionRHS2(double[] values) {
        for(int i = 0; i < this.convectionRHS2.length; ++i) {
            this.convectionRHS2[i] = values[i];
        }

    }

    public void setJointRHS1(double value) {
        this.jointRHS1 = value;
    }

    public void setJointRHS2(double value) {
        this.jointRHS2 = value;
    }

    private int setBiomass(double[] values, boolean delta) {
        int numLost = 0;

        for(int i = 0; i < this.biomass.length; ++i) {
            if (delta) {
                double[] var10000 = this.biomass;
                var10000[i] += values[i];
            } else if (values[i] < 0.0) {
                this.biomass[i] = 0.0;
            } else {
                this.biomass[i] = values[i];
            }

            if (this.biomass[i] < this.cParams.getMinSpaceBiomass()) {
                this.biomass[i] = 0.0;
                ++numLost;
            }
        }

        if (numLost == this.biomass.length) {
            return this.die();
        } else {
            if (this.cParams.showGraphics()) {
                this.cellColor = this.calculateColor();
            }

            this.updateDiffusibility();
            return 1;
        }
    }

    private int setBiomass3D(double[] values, boolean delta) {
        int numLost = 0;

        for(int i = 0; i < this.biomass.length; ++i) {
            if (delta) {
                double[] var10000 = this.biomass;
                var10000[i] += values[i];
            } else if (values[i] < 0.0) {
                this.biomass[i] = 0.0;
            } else {
                this.biomass[i] = values[i];
            }

            if (this.biomass[i] < this.cParams.getMinSpaceBiomass()) {
                this.biomass[i] = 0.0;
                ++numLost;
            }
        }

        if (numLost == this.biomass.length) {
            return this.die();
        } else {
            if (this.cParams.showGraphics()) {
                this.cellColor = this.calculateColor();
            }

            this.updateDiffusibility3D();
            return 1;
        }
    }

    private void updateDiffusibility() {
        int i;
        for(i = 0; i < this.biomass.length; ++i) {
            this.world.setDiffuseBiomassIn(this.x, this.y, i, true);
            if (this.biomass[i] == 0.0 && !this.cParams.allowCellOverlap()) {
                this.world.setDiffuseBiomassIn(this.x, this.y, i, false);
                this.world.setDiffuseBiomassOut(this.x, this.y, i, false);
            }
        }

        if (Utility.sum(this.biomass) >= this.cParams.getMaxSpaceBiomass()) {
            for(i = 0; i < this.biomass.length; ++i) {
                this.world.setDiffuseBiomassIn(this.x, this.y, i, false);
            }
        }

    }

    private void updateDiffusibility3D() {
        int i;
        for(i = 0; i < this.biomass.length; ++i) {
            this.world3D.setDiffuseBiomassIn(this.x, this.y, this.z, i, true);
            if (this.biomass[i] == 0.0 && !this.cParams.allowCellOverlap()) {
                this.world3D.setDiffuseBiomassIn(this.x, this.y, this.z, i, false);
                this.world3D.setDiffuseBiomassOut(this.x, this.y, this.z, i, false);
            }
        }

        if (Utility.sum(this.biomass) >= this.cParams.getMaxSpaceBiomass()) {
            for(i = 0; i < this.biomass.length; ++i) {
                this.world3D.setDiffuseBiomassIn(this.x, this.y, this.z, i, false);
            }
        }

    }

    public String attributesString() {
        return "temp attributes string";
    }

    private int calculateColor() {
        if (this.biomass.length == 2) {
            return Utility.pColor((int)this.biomass[1] * 25, (int)this.biomass[0] * 255 / 10, 55);
        } else {
            return this.biomass.length == 3 ? Utility.pColor((int)this.biomass[1] * 25, (int)this.biomass[0] * 255 / 10, 55 + (int)this.biomass[2] * 255 / 10) : Utility.pColor(0, (int)this.biomass[0] * 25, 55);
        }
    }

    public String cellType() {
        return "FBACell";
    }

    public int die() {
        return -1;
    }

    public int getColor() {
        return this.cellColor;
    }

    public synchronized double[] getBiomass() {
        return this.biomass;
    }

    public synchronized String[] getCellModelIDs() {
        String[] cellModelIDs = new String[this.fbaModels.length];

        for(int i = 0; i < this.fbaModels.length; ++i) {
            cellModelIDs[i] = this.fbaModels[i].getModelID();
        }

        return cellModelIDs;
    }

    public synchronized double[] getConvectionRHS1() {
        return this.convectionRHS1;
    }

    public synchronized double[] getConvectionRHS2() {
        return this.convectionRHS2;
    }

    public synchronized double getJointRHS2() {
        return this.jointRHS2;
    }

    public synchronized double getJointRHS1() {
        return this.jointRHS1;
    }

    public double[] getDeltaBiomass() {
        return this.deltaBiomass;
    }

    public int[] getFBAstatus() {
        return this.FBAstatus;
    }

    public boolean getStationaryStatus() {
        return this.stationaryStatus;
    }

    public void setStationaryStatus() {
        this.stationaryStatus = false;
    }

    public double[][] getFluxes() {
        return this.fluxes;
    }

    public int getID() {
        return this.id;
    }

    public int run() {
        return this.run(this.fbaModels);
    }

    public synchronized int run(FBAModel[] models) {
        this.deltaBiomass = new double[models.length];
        this.dyingBiomass = new double[this.biomass.length];
        this.allModelsGrowthRates = new double[models.length];
        this.deltaMedia = new double[models.length][];
        this.FBAstatus = new int[models.length];
        double biomassGrowthRate = 0.0;
        double rho = 1.0;
        double[][] allExchFluxes = new double[models.length][];
        double[][] lb = new double[models.length][];
        double[][] ub = new double[models.length][];

        int a;
        int i;
        double[] media;
        double[] rates;
        double[] var10000;
        double vMax;
        double[] mediaDelta;
        double maintenanceFlux;
        int i;
        double bioub;
        double km;
        for(a = 0; a < models.length; ++a) {
            i = a;
            if (this.biomass[a] != 0.0 && !(Utility.sum(this.biomass) >= this.cParams.getMaxSpaceBiomass())) {
                if ((!this.cParams.getSimulateActivation() || models[a].activate(this.cParams.getActivateRate())) && !this.stationaryStatus) {
                    media = null;
                    if (this.cParams.getNumLayers() == 1) {
                        media = this.world.getModelMediaAt(this.x, this.y, a);
                    } else if (this.cParams.getNumLayers() > 1) {
                        media = this.world3D.getModelMediaAt(this.x, this.y, this.z, a);
                    }

                    int j;
                    double maintenanceFlux;
                    lb[a] = models[a].getBaseExchLowerBounds();
                    ub[a] = models[a].getBaseExchUpperBounds();
                    this.applySignals(models[a], media);
                    rates = new double[lb[a].length];
                    label344:
                    switch (this.pParams.getExchangeStyle()) {
                        case MONOD:
                            double[] kmArr = models[a].getExchangeKm();
                            double[] vMaxArr = models[a].getExchangeVmax();
                            double[] hillCoeffArr = models[a].getExchangeHillCoefficients();
                            int j = 0;

                            while(true) {
                                if (j >= lb[i].length) {
                                    break label344;
                                }

                                if (lb[i][j] > 0.0) {
                                    rates[j] = -lb[i][j];
                                } else {
                                    km = FBAParameters.getDefaultKm();
                                    if (kmArr != null && kmArr.length > j && kmArr[j] > 0.0) {
                                        km = kmArr[j];
                                    }

                                    vMax = FBAParameters.getDefaultVmax();
                                    if (vMaxArr != null && vMaxArr.length > j && vMaxArr[j] > 0.0) {
                                        vMax = vMaxArr[j];
                                    }

                                    maintenanceFlux = FBAParameters.getDefaultHill();
                                    if (hillCoeffArr != null && hillCoeffArr.length > j && hillCoeffArr[j] > 0.0) {
                                        maintenanceFlux = hillCoeffArr[j];
                                    }

                                    if (media[j] / (this.cParams.getTimeStep() * this.biomass[i]) < this.calcMichaelisMentenRate(media[j] / this.cParams.getSpaceVolume(), km, vMax, maintenanceFlux)) {
                                        rates[j] = Math.min(Math.abs(lb[i][j]), Math.abs(media[j] / (this.cParams.getTimeStep() * this.biomass[i])));
                                    } else {
                                        rates[j] = Math.min(Math.abs(lb[i][j]), Math.abs(this.calcMichaelisMentenRate(media[j] / this.cParams.getSpaceVolume(), km, vMax, maintenanceFlux)));
                                    }
                                }

                                ++j;
                            }
                        case PSEUDO_MONOD:
                            double[] alphaArr = models[a].getExchangeAlphaCoefficients();
                            double[] wArr = models[a].getExchangeWCoefficients();
                            j = 0;

                            while(true) {
                                if (j >= lb[i].length) {
                                    break label344;
                                }

                                if (lb[i][j] > 0.0) {
                                    rates[j] = -lb[i][j];
                                } else {
                                    vMax = FBAParameters.getDefaultAlpha();
                                    if (alphaArr != null && alphaArr.length > j && alphaArr[j] > 0.0) {
                                        vMax = alphaArr[j];
                                    }

                                    maintenanceFlux = FBAParameters.getDefaultW();
                                    if (wArr != null && wArr.length > j && wArr[j] > 0.0) {
                                        maintenanceFlux = wArr[j];
                                    }

                                    rates[j] = Math.min(Math.abs(lb[i][j]), Math.abs(this.calcPseudoMonodRate(media[j] / this.cParams.getSpaceVolume(), vMax, maintenanceFlux)));
                                }

                                ++j;
                            }
                        default:
                            for(j = 0; j < lb[i].length; ++j) {
                                if (lb[i][j] > 0.0) {
                                    rates[j] = -lb[i][j];
                                } else {
                                    rates[j] = Math.min(Math.abs(lb[i][j]), Math.abs(this.calcStandardExchange(media[j])));
                                }
                            }
                    }

                    double[][] lightAbsorption = models[i].getLightAbsorption();

                    for(i = 0; i < lb[i].length; ++i) {
                        if (lightAbsorption[i][0] + lightAbsorption[i][1] > 0.0) {
                            rates[i] = Math.min(Math.abs(lb[i][i]), this.calcMaxLightUptake(media[i], this.biomass[i], this.cParams.getSpaceWidth(), lightAbsorption[i], this.cParams.getSpaceVolume()));
                        }
                    }

                    for(i = 0; i < lb[i].length; ++i) {
                        lb[i][i] = -1.0 * rates[i] / rho;
                    }

                    models[i].setExchLowerBounds(lb[i]);
                    bioub = (this.cParams.getMaxSpaceBiomass() - (Utility.sum(this.biomass) + Utility.sum(this.deltaBiomass))) / (this.biomass[i] * this.cParams.getTimeStep());
                    double currentbioub = models[i].getUpperBounds()[models[i].getBiomassReaction() - 1];
                    models[i].setBiomassUpperBound(Math.min(currentbioub, bioub));
                    j = models[i].run();
                    this.fluxes[i] = models[i].getFluxes();
                    double[] exchFlux = models[i].getExchangeFluxes();
                    allExchFluxes[i] = exchFlux;
                    mediaDelta = new double[exchFlux.length];
                    if (j != 5 && j != 180) {
                        this.deltaBiomass[i] = 0.0;
                        Arrays.fill(mediaDelta, 0.0);
                        this.deltaMedia[i] = mediaDelta;
                    } else {
                        for(int j = 0; j < mediaDelta.length; ++j) {
                            if (lightAbsorption[j][0] + lightAbsorption[j][1] > 0.0) {
                                mediaDelta[j] = 0.0;
                            } else {
                                mediaDelta[j] = exchFlux[j] * this.biomass[i] * this.cParams.getTimeStep();
                            }
                        }

                        this.deltaMedia[i] = mediaDelta;
                        biomassGrowthRate = models[i].getBiomassFluxSolution();
                        this.deltaBiomass[i] = models[i].getBiomassFluxSolution() * this.cParams.getTimeStep() * this.biomass[i];
                        this.allModelsGrowthRates[i] = biomassGrowthRate;
                        var10000 = this.deltaBiomass;
                        var10000[i] *= 1.0 - models[i].getGenomeCost();
                        if (models[i].hasMaintenance()) {
                            maintenanceFlux = models[i].getMaintenanceFluxSolution();
                            if (maintenanceFlux < models[i].getMinMaintenanceFlux()) {
                                var10000 = this.deltaBiomass;
                                var10000[i] -= this.biomass[i] * (1.0 - maintenanceFlux / models[i].getMinMaintenanceFlux());
                            }
                        }

                        if (!this.pParams.getAllowFluxWithoutGrowth() && this.deltaBiomass[i] < 0.0) {
                            this.deltaBiomass[i] = 0.0;
                            Arrays.fill(this.deltaMedia[i], 0.0);
                        }

                        if (this.cParams.showGraphics()) {
                            this.cellColor = this.calculateColor();
                        }
                    }

                    Object[] temp = this.calcDeathRateAndMetConsumption(models[i], media, this.biomass[i]);
                    maintenanceFlux = (Double)temp[0];
                    Map<Integer, Double> consumed_mets = (Map)temp[1];
                    this.dyingBiomass[i] = maintenanceFlux;
                    Set<Integer> consumed_met_keys = consumed_mets.keySet();

                    int key;
                    for(Iterator var26 = consumed_met_keys.iterator(); var26.hasNext(); var10000[key] -= (Double)consumed_mets.get(key)) {
                        key = (Integer)var26.next();
                        var10000 = this.deltaMedia[i];
                    }
                }
            } else {
                this.deltaBiomass[a] = 0.0;
                this.dyingBiomass[a] = 0.0;
                media = models[a].getExchangeFluxes();
                rates = new double[media.length];
                Arrays.fill(rates, 0.0);
                this.deltaMedia[a] = rates;
            }
        }

        if (!this.stationaryStatus) {
            double[][] uptakeMat = this.world.simulateCellUpdateMedia(this.x, this.y, models, this.deltaMedia);
            boolean reOptimizeFlag = false;
            media = this.world.getMediaAt(this.x, this.y);
            rates = new double[media.length];

            int a;
            int j;
            for(a = 0; a < media.length; ++a) {
                bioub = 0.0;
                ArrayList<Integer> uptakingModels = new ArrayList();

                for(int l = 0; l < models.length; ++l) {
                    bioub += uptakeMat[l][a];
                    if (uptakeMat[l][a] < 0.0) {
                        uptakingModels.add(l);
                    }
                }

                if (bioub < 0.0 && bioub < -media[a]) {
                    reOptimizeFlag = true;

                    Integer l;
                    for(Iterator var40 = uptakingModels.iterator(); var40.hasNext(); lb[l][j] = -vMax / (this.biomass[l] * this.cParams.getTimeStep())) {
                        l = (Integer)var40.next();
                        vMax = media[a] * (uptakeMat[l][a] / bioub);
                        int[] modelMediaIndexes = this.world.getModelMediaIndexes(this.x, this.y, l);
                        j = ArrayUtils.indexOf(modelMediaIndexes, a);
                    }
                }
            }

            if (reOptimizeFlag) {
                for(a = 0; a < models.length; ++a) {
                    i = a;
                    if (this.biomass[a] != 0.0 && !(Utility.sum(this.biomass) >= this.cParams.getMaxSpaceBiomass()) && (!this.cParams.getSimulateActivation() || models[a].activate(this.cParams.getActivateRate())) && !this.stationaryStatus) {
                        models[a].setExchLowerBounds(lb[a]);
                        double bioub = (this.cParams.getMaxSpaceBiomass() - (Utility.sum(this.biomass) + Utility.sum(this.deltaBiomass))) / (this.biomass[a] * this.cParams.getTimeStep());
                        km = models[a].getUpperBounds()[models[a].getBiomassReaction() - 1];
                        models[a].setBiomassUpperBound(Math.min(km, bioub));
                        int stat = models[a].run();
                        this.fluxes[a] = models[a].getFluxes();
                        if (stat != 5 && stat != 180) {
                            this.deltaBiomass[a] = 0.0;
                        }

                        if (stat == 5 || stat == 180) {
                            mediaDelta = models[a].getExchangeFluxes();
                            allExchFluxes[a] = mediaDelta;
                            double[] mediaDelta = new double[mediaDelta.length];

                            for(j = 0; j < mediaDelta.length; ++j) {
                                mediaDelta[j] = mediaDelta[j] * this.biomass[i] * this.cParams.getTimeStep();
                                double[][] lightAbsorption = models[i].getLightAbsorption();
                                if (lightAbsorption[j][0] + lightAbsorption[j][1] > 0.0) {
                                    mediaDelta[j] = 0.0;
                                } else {
                                    mediaDelta[j] = mediaDelta[j] * this.biomass[i] * this.cParams.getTimeStep();
                                }
                            }

                            this.deltaMedia[i] = mediaDelta;
                            biomassGrowthRate = models[i].getBiomassFluxSolution();
                            this.deltaBiomass[i] = models[i].getBiomassFluxSolution() * this.cParams.getTimeStep() * this.biomass[i];
                            this.allModelsGrowthRates[i] = biomassGrowthRate;
                            var10000 = this.deltaBiomass;
                            var10000[i] *= 1.0 - models[i].getGenomeCost();
                            if (models[i].hasMaintenance()) {
                                maintenanceFlux = models[i].getMaintenanceFluxSolution();
                                if (maintenanceFlux < models[i].getMinMaintenanceFlux()) {
                                    var10000 = this.deltaBiomass;
                                    var10000[i] -= this.biomass[i] * (1.0 - maintenanceFlux / models[i].getMinMaintenanceFlux());
                                }
                            }

                            if (!this.pParams.getAllowFluxWithoutGrowth() && this.deltaBiomass[i] < 0.0) {
                                this.deltaBiomass[i] = 0.0;

                                for(j = 0; j < this.deltaMedia[i].length; ++j) {
                                    this.deltaMedia[i][j] = 0.0;
                                }
                            }
                        }
                    }
                }
            }

            for(a = 0; a < models.length; ++a) {
                if (this.biomass[a] != 0.0 && !(Utility.sum(this.biomass) >= this.cParams.getMaxSpaceBiomass())) {
                    if (this.cParams.getNumLayers() == 1) {
                        this.world.changeModelMedia(this.x, this.y, a, this.deltaMedia[a]);
                    } else if (this.cParams.getNumLayers() > 1) {
                        this.world3D.changeModelMedia(this.x, this.y, this.z, a, this.deltaMedia[a]);
                    }
                }
            }
        }

        if (this.cParams.getBatchDilution()) {
            this.stationaryStatus = true;

            for(a = 0; a < this.deltaMedia.length; ++a) {
                if (this.deltaMedia[a] != null) {
                    for(i = 0; i < this.deltaMedia[a].length; ++i) {
                        if (this.deltaMedia[a][i] != 0.0) {
                            this.stationaryStatus = false;
                        }
                    }
                }
            }
        }

        if (this.cParams.showGraphics()) {
            this.cellColor = this.calculateColor();
        }

        System.out.println("TESTING 3: " + Arrays.toString(this.deltaBiomass));
        return this.updateCellData(this.deltaBiomass, this.fluxes, this.allModelsGrowthRates);
    }

    private double calcMichaelisMentenRate(double mediaConc, double km, double vMax, double hill) {
        return mediaConc * vMax / (km + mediaConc);
    }

    private double calcPseudoMonodRate(double mediaConc, double alpha, double w) {
        return Math.min(mediaConc * alpha, w);
    }

    private double calcStandardExchange(double mediaConc) {
        return mediaConc;
    }

    private double addDemographicNoise(double currentBiomass, double biomassGrowthRate, double demographicNoiseSigmaZero) {
        double noisyBiomass = currentBiomass;
        if (biomassGrowthRate > 0.0) {
            double noiseSigma = biomassGrowthRate * demographicNoiseSigmaZero;
            double poissonLambda = 2.0 * currentBiomass / (this.cParams.getTimeStep() * noiseSigma * noiseSigma);
            if (poissonLambda > 0.0) {
                double gammaAlpha = Poisson.random(poissonLambda, new MersenneTwister());
                if (gammaAlpha > 0.0) {
                    double gammaSample = Gamma.random(gammaAlpha, 1.0, new MersenneTwister());
                    noisyBiomass = 0.5 * gammaSample * this.cParams.getTimeStep() * noiseSigma * noiseSigma;
                } else if (gammaAlpha == 0.0) {
                    noisyBiomass = 0.0;
                }
            } else if (poissonLambda == 0.0) {
                noisyBiomass = 0.0;
            }
        }

        return noisyBiomass;
    }

    private Object[] calcDeathRateAndMetConsumption(FBAModel model, double[] media, double biomass) {
        double death_rate = 0.0;
        Map<Integer, Double> consumed_mets = new HashMap();
        double space_volume = this.cParams.getSpaceVolume();
        Iterator var10 = model.getSignals().iterator();

        while(var10.hasNext()) {
            Signal signal = (Signal)var10.next();
            if (signal.affectsDeathRate()) {
                int signal_met = signal.getExchMet();
                double death_caused_by_toxin = signal.calculateDeathRate(media[signal_met] / space_volume);
                death_caused_by_toxin = death_caused_by_toxin * biomass * this.cParams.getTimeStep();
                death_rate += death_caused_by_toxin;
                if (signal.isMetConsumed()) {
                    consumed_mets.put(signal_met, death_caused_by_toxin);
                }
            }
        }

        return new Object[]{death_rate, consumed_mets};
    }

    private boolean applySignals(FBAModel model, double[] media) {
        if (model.getSignals().size() > 0) {
            double[] all_lb = model.getLowerBounds();
            double[] all_ub = model.getUpperBounds();
            double space_volume = this.cParams.getSpaceVolume();
            Iterator var7 = model.getSignals().iterator();

            while(var7.hasNext()) {
                Signal signal = (Signal)var7.next();
                if (signal.getReaction() != -1) {
                    int signal_met;
                    int signal_rxn;
                    double new_ub;
                    if (signal.affectsLb()) {
                        signal_met = signal.getExchMet();
                        signal_rxn = signal.getReaction();
                        new_ub = signal.calculateBound(media[signal_met] / space_volume);
                        all_lb[signal_rxn] = new_ub;
                    }

                    if (signal.affectsUb()) {
                        signal_met = signal.getExchMet();
                        signal_rxn = signal.getReaction();
                        new_ub = signal.calculateBound(media[signal_met] / space_volume);
                        all_ub[signal_rxn] = new_ub;
                    }
                }
            }

            model.setLowerBounds(all_lb);
            model.setUpperBounds(all_ub);
        }

        return true;
    }

    private double calcMaxLightUptake(double lightFlux, double biomass, double gridSize, double[] absorption, double gridVolume) {
        double biomassAbsorption = absorption[1] * biomass / (gridVolume * 1.0E-6);
        double absorbance = (absorption[0] + biomassAbsorption) * gridSize * 0.01;
        double deltaFlux = lightFlux * (1.0 - Math.exp(-absorbance));
        double absorbedFlux = deltaFlux * biomassAbsorption / (biomassAbsorption + absorption[0]);
        double absorbedPhotonsPerHourPerBiomass = 0.36000000000000004 * absorbedFlux * gridSize * gridSize / biomass;
        return absorbedPhotonsPerHourPerBiomass;
    }

    public int updateCellData(double[] deltaBiomass, double[][] fluxes, double[] biomassGrowthRates) {
        this.deltaBiomass = deltaBiomass;
        System.out.println("TESTING 2: " + Arrays.toString(deltaBiomass));
        this.fluxes = fluxes;
        int numDead = 0;

        for(int i = 0; i < this.biomass.length; ++i) {
            double[] var10000 = this.dyingBiomass;
            var10000[i] += this.cParams.getDeathRate() * this.biomass[i] * this.cParams.getTimeStep();
            var10000 = this.biomass;
            var10000[i] += deltaBiomass[i];
            var10000 = this.biomass;
            var10000[i] -= this.dyingBiomass[i];
            if (this.fbaModels[i].getNeutralDrift() && deltaBiomass[i] > 0.0 && this.cParams.getDeathRate() == 0.0) {
                this.biomass[i] = this.addDemographicNoise(this.biomass[i], biomassGrowthRates[i], this.fbaModels[i].getNeutralDriftSigma());
            } else if (this.fbaModels[i].getNeutralDrift() && this.cParams.getDeathRate() != 0.0) {
                System.out.println("Error in model " + i + ": Demographic noise is applies only if the death rate for the model is zero. Noise will not be applied.");
                System.err.println("Error in model " + i + ": Demographic noise is applies only if the death rate for the model is zero. Noise will not be applied.");
            }

            if (this.biomass[i] < this.cParams.getMinSpaceBiomass()) {
                this.biomass[i] = 0.0;
            }

            if (this.biomass[i] == 0.0) {
                ++numDead;
            }
        }

        if (this.cParams.showGraphics()) {
            this.cellColor = this.calculateColor();
        }

        if (numDead == this.biomass.length) {
            return this.die();
        } else {
            if (this.cParams.getNumLayers() == 1) {
                this.updateDiffusibility();
            } else if (this.cParams.getNumLayers() > 1) {
                this.updateDiffusibility3D();
            }

            return 1;
        }
    }

    public String statisticsString() {
        return "temp statistics string";
    }

    public String toString() {
        String str = "FBA Cell\nPos: (" + this.x + ", " + this.y + ")\nBiomass: " + this.biomass[0];

        int i;
        for(i = 1; i < this.biomass.length; ++i) {
            str = str + ", " + this.biomass[i];
        }

        str = str + "\n";
        if (this.fluxes != null) {
            for(i = 0; i < this.fluxes.length; ++i) {
                str = str + "model " + i + ": ";
                if (this.fluxes[i] == null) {
                    str = str + "no current fluxes";
                } else {
                    str = str + "fluxes: ";

                    for(int j = 0; j < this.fluxes[i].length; ++j) {
                        str = str + this.fluxes[i][j] + "\t";
                    }
                }

                str = str + "\n";
            }
        }

        return str;
    }

    public Cell backup(World2D backupWorld) {
        FBACell bak = new FBACell(this.x, this.y, (double[])this.biomass.clone(), (FBAWorld)backupWorld, this.fbaModels, this.cParams, this.pParams);
        return bak;
    }

    public Cell backup(World3D backupWorld) {
        FBACell bak = new FBACell(this.x, this.y, this.z, (double[])this.biomass.clone(), (FBAWorld3D)backupWorld, this.fbaModels, this.cParams, this.pParams);
        return bak;
    }

    public void changeModelsInCell(Model[] oldModels, Model[] newModels) {
        int[] idx = new int[oldModels.length];

        int i;
        for(i = 0; i < oldModels.length; ++i) {
            idx[i] = -1;
        }

        for(i = 0; i < oldModels.length; ++i) {
            for(int j = 0; j < newModels.length; ++j) {
                if (oldModels[i].equals(newModels[j])) {
                    idx[i] = j;
                }
            }
        }

        double[] newBiomass = new double[newModels.length];
        double[] newDeltaBiomass = new double[newModels.length];
        double[][] newFluxes = new double[newModels.length][];

        int i;
        for(i = 0; i < idx.length; ++i) {
            if (idx[i] != -1) {
                newBiomass[idx[i]] = this.biomass[i];
                newDeltaBiomass[idx[i]] = this.deltaBiomass[i];
                newFluxes[idx[i]] = this.fluxes[i];
            }
        }

        this.biomass = newBiomass;
        this.deltaBiomass = newDeltaBiomass;
        this.fluxes = newFluxes;
        this.fbaModels = new FBAModel[newModels.length];

        for(i = 0; i < newModels.length; ++i) {
            this.fbaModels[i] = (FBAModel)newModels[i];
        }

    }

    public CometsParameters getCometsParameters() {
        return this.cParams;
    }

    public double[] getMedia() {
        if (this.cParams.getNumLayers() == 1) {
            return this.world.getMediaAt(this.x, this.y);
        } else {
            return this.cParams.getNumLayers() > 1 ? this.world3D.getMediaAt(this.x, this.y, this.z) : null;
        }
    }

    public Comets getComets() {
        if (this.cParams.getNumLayers() == 1) {
            return this.world.getComets();
        } else {
            return this.cParams.getNumLayers() > 1 ? this.world3D.getComets() : null;
        }
    }

    public void refreshParameters() {
        Comets c = null;
        if (this.world != null) {
            c = this.world.getComets();
        } else {
            c = this.world3D.getComets();
        }

        this.pParams = (FBAParameters)c.getPackageParameters();
        this.cParams = c.getParameters();
    }

    public void setParameters(CometsParameters cParams) {
        this.cParams = cParams;
    }

    public void setPackageParameters(PackageParameters pParams) {
        if (pParams.getClass().equals(FBAParameters.class)) {
            this.pParams = (FBAParameters)pParams;
        }

    }
}
