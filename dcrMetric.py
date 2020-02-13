import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats


class DcrMetric:
    """Parameterize the observing conditions into a single `visit_measure`,
    and evaluate how well DCR is constrained compared to a reference
    set of observations.

    Attributes
    ----------
    density : `scipy.stats.kde.gaussian_kde`
        Raw KDE of the visit measures, used as the weights for `kde`.
    kde : `scipy.stats.kde.gaussian_kde`
        Normalized Kernel Density Estimate (KDE) of the visit measures.
    max_airmass : `float`
        The maximum airmass to include when evaluating the KDE.
    metric : `float`
        The metric evaluating how well DCR is constrained for the given KDE.
    metric_reference : `float`
        The `metric` calculated for `ref_kde`
    norm : `float`
        The normalization factor required for
        the KDE to integrate to 1.0 over `range`.
    range : `tuple` of two `float`s.
        The minimum and maximum visit measures, derived from `max_airmass`.
    ref_density_kde : `scipy.stats.kde.gaussian_kde`
        Raw KDE of the reference set of visit measures.
        Used as the weights for `ref_kde`.
    ref_kde : `scipy.stats.kde.gaussian_kde`
        Normalized KDE of the reference set of visit measures.
    dataframe : `pandas.DataFrame`
        The record of the airmass and hour angle of each input observation.
    width : `float`
        Width of the gaussians used as the basis for the KDEs.

    Notes
    -----
    The metric is calculated as the significance of each input observation,
    added in quadrature. The significance of an observation is the the inverse
    of the weight, or density, KDE evaluated at the observation's visit measure.
    Thus, the significance of an individual observation is 1 when it is the
    only one constraining its region, but decreases when there are many other
    observations at similar observing conditions.

    For ease of use, the final metric is divided by the metric calculated for
    a set of well-sampled reference observations, so that any metric greater
    than one signifies that the model is sufficiently constrained.
    """

    def __init__(self, max_airmass=1.8, width=0.4):
        self.max_airmass = max_airmass
        self.range = (1 - max_airmass, max_airmass - 1)
        self.width = width
        self.kde = None
        self.density = None
        self.norm = None
        self.metric = None
        # visit will be the index of the dataframe.
        self.dataframe = pd.DataFrame(columns=['airmass', 'hour_angle'], dtype='float')
        self.metric_reference = self.calculateReference()

    def parameterizeVisit(self, airmass, hour_angle):
        """Convert airmass and hour angle to visit measure.

        Parameters
        ----------
        airmass: `float`
            The observed airmass of the observation.
        hour_angle: `float`
            The observed hour angle of the observation.

        Returns
        -------
        visit_measure: `float`
            A one dimensional parameterization of airmass and hour angle.
        """
        visit_measure = (airmass - 1.)*np.sign(hour_angle)
        return visit_measure

    def parameterizeDataframe(self, visits=None):
        """Calculate the visit measures for all observations in the dataframe.

        Parameters
        ----------
        visits: `list`, optional
            The visit identifiers of the desired observations.
            If `None`, all observations in the dataframe are used.

        Returns
        -------
        visit_measures : `list` of `float`
            The one-dimensional parameterization of airmass and hour angle
            for all visits in the dataframe.
        """
        visit_measures = []
        if visits is None:
            visits = self.dataframe.index
        for visit in visits:
            airmass = self.dataframe['airmass'][visit]
            hour_angle = self.dataframe['hour_angle'][visit]
            measure = self.parameterizeVisit(airmass, hour_angle)
            if measure >= self.range[0] and measure <= self.range[1]:
                visit_measures.append(measure)
        return visit_measures

    def addVisit(self, visit, airmass, hour_angle):
        """Add a new visit to the dataframe.

        Parameters
        ----------
        visit: `float`
            The visit identifier of the observation to be added.
        airmass: `float`
            The observed airmass of the observation.
        hour_angle: `float`
            The observed hour angle of the observation.
        """
        if visit in self.dataframe.index:
            raise KeyError(f"Visit {visit} is already present in the dataframe,"
                           " and cannot be added. Skipping it.")
        self.dataframe.loc[visit] = [airmass, hour_angle]
        self.kde = None

    def addVisitInfo(self, visitInfo):
        """Add a new visit to the dataframe using the VisitInfo of an exposure.

        Parameters
        ----------
        visitInfo: `lsst.afw.image.VisitInfo`
            The VisitInfo metadata of an exposure.
        """
        visit = visitInfo.getExposureId()//100  # This is Decam-specific
        hour_angle = visitInfo.getBoresightHourAngle().asRadians()
        airmass = visitInfo.getBoresightAirmass()
        self.addVisit(visit, airmass, hour_angle)

    def delVisit(self, visit):
        """Remove a visit from the dataframe.

        Parameters
        ----------
        visit: `float`
            The visit identifier of the observation to be removed.
        """
        try:
            self.dataframe.drop(visit, inplace=True)
        except Exception:
            raise KeyError(f"Could not remove visit {visit}, it is not in the dataframe.")
        else:
            self.kde = None

    def testVisit(self, airmass=None, hour_angle=None, visit_measure=None):
        """Evaluate how significant a new observation will be to constrain DCR.

        Parameters
        ----------
        airmass: `float`, optional
            The observed airmass of the observation.
        hour_angle: `float`, optional
            The observed hour angle of the observation.
        visit_measure: `float`, optional
            The one-dimensional parameterization of airmass and hour angle.
        """
        self.checkKde()
        if visit_measure is None:
            visit_measure = self.parameterizeVisit(airmass, hour_angle)
        weight = 1./(self.evaluateKde(visit_measure)*self.density(visit_measure))
        return weight

    def checkKde(self):
        """Check if the KDE exists, and calculate it from the dataframe if not.
        """
        if self.kde is None:
            self.kde, self.density = self.calculateKde(return_weights=True)
            self.norm = self.calculateNormalization()
            self.metric = self.calculateMetric()

    def calculateKde(self, visit_measures=None, return_weights=False):
        """Calculate the KDE for a given list of visit measures.

        Parameters
        ----------
        visit_measures: `list` of `float`, optional
            The list of visit measures to calculate the KDE for.
            If not supplied, the visit measures stored
            in the dataframe will be used.
        return_weights: `bool`
            Return the weights KDE as well as the normalized KDE?

        Returns
        -------
        kde: `scipy.stats.kde.gaussian_kde`
            Normalized KDE of the visit measures.
        density_kde: `scipy.stats.kde.gaussian_kde`
            Raw KDE of the visit measures.
            Only returned if ``return_weights`` is set.
        """
        if visit_measures is None:
            visit_measures = self.parameterizeDataframe()
        density_kde = scipy.stats.gaussian_kde(visit_measures, bw_method=self.width**2)
        visit_weights = 1/density_kde(visit_measures)
        kde = scipy.stats.gaussian_kde(visit_measures, bw_method=self.width, weights=visit_weights)
        if return_weights:
            return(kde, density_kde)
        else:
            return kde

    def calculateReference(self):
        """Calculate the metric for a reference distribution of visit measures.

        Returns
        -------
        metric: `float`
            The metric of the reference visit measures.
        """
        test_measures = np.linspace(self.range[0], self.range[1], num=8, endpoint=True)
        kde, density_kde = self.calculateKde(visit_measures=test_measures, return_weights=True)
        self.ref_kde = kde
        self.ref_density_kde = density_kde
        metric = self._calculateMetric(density_kde=density_kde, visit_measures=test_measures,
                                       use_reference=False)
        return metric

    def calculateNormalization(self, kde=None):
        """Normalize the KDE to integrate to one over one unit of visit measure.

        Parameters
        ----------
        kde: `scipy.stats.kde.gaussian_kde`, optional
            The KDE to normalize.
            Uses ``self.kde`` if not supplied.

            Returns
        -------
        normalization: `float`
            The normalization factor for the given KDE.
        """
        if kde is None:
            self.checkKde()
            kde = self.kde
        normalization = (self.range[1] - self.range[0])/kde.integrate_box(*self.range)
        return normalization

    def calculateMetric(self, kde=None):
        """Calculate the DCR metric for a given KDE.

        Parameters
        ----------
        kde: `scipy.stats.kde.gaussian_kde`, optional
            Uses ``self.kde`` if not supplied.

        Returns
        -------
        metric: `float`
            The metric evaluating how well DCR is constrained for the given KDE.
        """
        if kde is None:
            self.checkKde()
            kde = self.kde
        metric = self._calculateMetric(self.density, self.kde.dataset)
        return metric

    def _calculateMetric(self, density_kde=None, visit_measures=None, use_reference=True):
        if density_kde is None:
            density_kde = self.density
        if visit_measures is None:
            visit_measures = self.kde.dataset
        visit_weights = 1./(density_kde(visit_measures))**2
        metric = np.sqrt(np.sum(visit_weights))
        if use_reference:
            metric /= self.metric_reference
        return metric

    def evaluateKde(self, visit_measure):
        """Evaluate how well DCR is constrained for an observation.

        Parameters
        ----------
        visit_measure: `float`
            The one-dimensional parameterization of airmass and hour angle.

        Returns
        -------
        visit_metric: `float`
            The metric evaluating how well a new observation
            can be constrained for the current DCR model.
        """
        template_metric = self.metric
        visit_metric = self.kde(visit_measure)
        visit_metric *= self.norm*template_metric
        return visit_metric

    def evaluateVisit(self, airmass, hour_angle):
        """Evaluate how well DCR is constrained for an observation.

        Parameters
        ----------
        airmass: `float`
            The observed airmass of the observation.
        hour_angle: `float`
            The observed hour angle of the observation.

        Returns
        -------
        visit_metric: `float`
            The metric evaluating how well a new observation
            can be constrained for the current DCR model.
        """
        self.checkKde()
        visit_measure = self.parameterizeVisit(airmass, hour_angle)
        visit_metric = self.evaluateKde(visit_measure)
        return visit_metric

    def evaluateVisitInfo(self, visitInfo):
        """Evaluate how well DCR is constrained for an observation.

        Parameters
        ----------
        visitInfo: `lsst.afw.image.VisitInfo`
            The VisitInfo metadata of an exposure.

        Returns
        -------
        visit_metric: `float`
            The metric evaluating how well a new observation
            can be constrained for the current DCR model.
        """
        hour_angle = visitInfo.getBoresightHourAngle().asRadians()
        airmass = visitInfo.getBoresightAirmass()
        return self.evaluateVisit(airmass, hour_angle)

    def plotMetric(self, window=1000, plot_input_observations=True):
        """Generate a plot to visualize the KDE and the input visit measures.

        Parameters
        ----------
        window: `int`
            The figure number to add the new plots to.
        plot_input_observations: `bool`
            Overplot the visit measures of the input observations.

        Notes
        -----
        Visit measure is a parameterization of airmass and hour angle, which
        takes advantage of the fact that, for a single observatory, a given
        target will trace the same arc across the sky. Thus, the observing
        conditions can be uniquely defined by the airmass and whether the
        target is rising or setting, i.e. the sign of the hour angle.
        Since the minimum airmass is 1.0, we define the visit measure as
        (airmass - 1) x sign(hour angle)
        """
        self.checkKde()
        visit_measures = self.parameterizeDataframe()
        xv = np.arange(self.range[0]-.3, self.range[1]+.3, .01)
        nVisits = len(visit_measures)
        plt.figure(window)
        plt.xlabel('Visit measure', size=12)
        plt.ylabel('DCR metric', size=12)
        metric_values = self.evaluateKde(xv)
        plt.plot(xv, metric_values, '-')
        if plot_input_observations:
            plt.plot(visit_measures, np.ones(nVisits), '+')
