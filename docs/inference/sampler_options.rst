.. _inference_sampler_options:

##################################
Choosing and configuring a sampler
##################################

``pycbc_inference`` can drive a number of different sampling engines. Which
one to use, and how to configure it, is set by the ``[sampler]`` section of
the configuration file. This page explains what the choices are and what
the options in that section do.

For ready-to-run configurations for each sampler, see
:ref:`inference_example_samplers`. For the internals, see the
:doc:`sampler API <sampler_api>`.

.. note::
   Most options in the ``[sampler]`` section are passed straight through to
   the underlying package, so its own documentation is the authoritative
   reference for what a given option means. The links in each section below
   point at it.

*******************
Which sampler?
*******************

The samplers fall into two families, and the family matters more than the
individual choice:

**Nested samplers** integrate the likelihood over the prior. They return an
estimate of the evidence (:math:`\log Z`) as well as posterior samples, they
decide for themselves when they have converged, and they need no burn-in or
autocorrelation analysis. This is usually what you want for a
gravitational-wave parameter estimation run.

**MCMC samplers** produce a chain of samples from the posterior. They do not
give you an evidence, and you have to tell them how long to run and discard
the burn-in yourself, but they can be more efficient for some problems and
are easier to reason about when something goes wrong.

.. list-table::
   :widths: 15 10 45
   :header-rows: 1

   * - Name
     - Type
     - Notes
   * - ``dynesty``
     - nested
     - The usual default. Well tested, checkpoints, many tuning options.
   * - ``ultranest``
     - nested
     - Robust when the posterior is multi-modal. Good diagnostics.
   * - ``nessai``
     - nested
     - Normalising flows; aimed at expensive likelihoods.
   * - ``multinest``
     - nested
     - Long established. Needs the external MultiNest library.
   * - ``cpnest``
     - nested
     - Parallel nested sampling.
   * - ``snowline``
     - importance
     - Cheap, but best on well behaved posteriors.
   * - ``emcee``
     - MCMC
     - Affine-invariant ensemble sampler.
   * - ``ptemcee``
     - MCMC
     - Parallel-tempered ``emcee``. Better when multi-modal.
   * - ``epsie``
     - MCMC
     - Parallel-tempered, with per-parameter jump proposals.

**********************
Options for every run
**********************

``name``
    Required, and must match one of the names above.

Two more controls are worth knowing about even though they are *not* set in
the ``[sampler]`` section:

* ``--nprocesses`` on the ``pycbc_inference`` command line sets the number of
  parallel processes. Every sampler on this page can use it.
* ``--seed`` sets the random seed. Note that not every underlying package
  seeds itself from this, so a fixed seed does not always give a bit-for-bit
  identical run.

*******************
Nested samplers
*******************

dynesty
=======

See the `dynesty documentation <https://dynesty.readthedocs.io/>`_.

``nlive``
    Required. The number of live points. More live points give a more finely
    sampled posterior and a more accurate evidence, at proportionally more
    cost. A **negative** value selects dynesty's dynamic nested sampler,
    which chooses the number of live points itself.

Stopping conditions
-------------------

``dlogz``
    Target remaining evidence. The usual way to control how long the run
    goes on; smaller runs longer.

``maxiter``, ``maxcall``
    Hard caps on iterations and likelihood calls.

``logl_max``
    Stop once the likelihood reaches this value.

``n_effective``
    Target effective sample size.

How points are sampled
----------------------

``sample``
    ``unif``, ``rwalk``, ``rwalk2`` (a modified random walk provided by
    pycbc), ``slice`` or ``rslice``. Uniform sampling is only viable in low
    dimensions; the walk and slice methods scale better.

``walks``
    For the random-walk methods, the minimum number of steps taken before a
    new point is accepted. Raising it decorrelates the chain at linear cost.

``maxmcmc``, ``nact``
    Used by the ``rwalk2`` method: the maximum number of steps taken when
    evolving a point, and the number of autocorrelation lengths required
    before stopping.

``bound``
    How the prior volume is bounded: ``none``, ``single``, ``multi``,
    ``balls`` or ``cubes``. ``multi`` handles multi-modal posteriors.

``update_interval``, ``first_update_min_ncall``, ``first_update_min_eff``
    Control how often the bounds are rebuilt, and how long the sampler waits
    before building the first one.

``enlarge``, ``bootstrap``
    Expand the bounding volume, either by a fixed factor or by
    bootstrapping. Note that dynesty does not accept both at once unless one
    is disabled.

Checkpointing
-------------

``checkpoint_time_interval``
    Seconds between checkpoints. If unset, the run does not checkpoint.

``maxcall``
    Also acts as how often the sampler checks whether to checkpoint.

ultranest
=========

See the `UltraNest documentation
<https://johannesbuchner.github.io/UltraNest/>`_.

``min_num_live_points``
    The number of live points.

``dlogz``, ``min_ess``, ``frac_remain``, ``max_iters``, ``max_ncalls``, ``max_num_improvement_loops``
    Stopping conditions.

``stepsampling``
    Turn on step sampling, which scales to higher dimension than the default
    region sampling.

``update_interval_iter_fraction``, ``update_interval_ncall``
    How often the bounding region is updated.

``log_dir``, ``log_interval``, ``show_status``, ``enable_plots``
    Diagnostic output.

nessai
======

See the `nessai documentation <https://nessai.readthedocs.io/>`_.

``nlive``
    The number of live points.

``importance_nested_sampler``
    Set to ``True`` to use the importance nested sampler (i-nessai) instead
    of the standard one. This changes which of the options below are valid,
    since the two samplers take different keywords.

Any keyword accepted by ``nessai``'s ``FlowSampler`` or its ``run`` method
may be given in this section; unrecognised options are rejected with an
error listing them.

snowline
========

See the `snowline documentation <https://johannesbuchner.github.io/snowline/>`_.

``num_global_samples``, ``num_gauss_samples``
    How many samples are drawn globally and from the fitted Gaussian.

``min_ess``
    Target effective sample size.

``max_ncalls``, ``max_improvement_loops``
    Stopping conditions.

multinest
=========

See the `PyMultiNest documentation
<https://johannesbuchner.github.io/PyMultiNest/>`_.

``nlivepoints``
    The number of live points.

``evidence-tolerance``
    Target accuracy on the evidence; the main stopping condition.

``sampling-efficiency``
    Trades sampling efficiency against accuracy. Lower values are more
    accurate and slower.

cpnest
======

See the `cpnest documentation <https://github.com/johnveitch/cpnest>`_.

``nlive``
    The number of live points.

``maxmcmc``
    The maximum length of the internal MCMC used to evolve a point.

``nthreads``
    Number of parallel threads cpnest should use.

*******************
MCMC samplers
*******************

Options shared by the MCMC samplers
===================================

``niterations``
    Run for a fixed number of iterations.

``effective-nsamples``
    Alternatively, run until this many *effective* samples have been
    collected, i.e. after burn-in and thinning by the autocorrelation
    length. Use this rather than ``niterations`` if you want the run to
    decide for itself when it has enough.

``max-samples-per-chain``
    Caps memory use by thinning on the fly once a chain reaches this length.
    Needed for long runs.

``checkpoint-interval``
    Iterations between checkpoints.

``checkpoint-signal``
    Signal to checkpoint and exit on, for running under a batch system.

``thin-interval``
    Thin the stored chain by this factor as it is written.

Burn-in is configured in its own ``[sampler-burn_in]`` section:

``burn-in-test``
    Which test decides that the chain has burned in, e.g. ``max_posterior``,
    ``halfchain``, ``ks_test``, or a combination such as
    ``nacl & max_posterior``.

emcee
=====

See the `emcee documentation <https://emcee.readthedocs.io/>`_.

``nwalkers``
    Required. The number of walkers in the ensemble. It must exceed twice
    the number of varied parameters, and more walkers explore better.

``logpost-function``
    Which model attribute to use as the log posterior.

ptemcee
=======

See the `ptemcee documentation <https://github.com/willvousden/ptemcee>`_.

``nwalkers``
    Required, as for ``emcee``.

``ntemps``, ``tmax``, ``betas``, ``betas-file``
    The temperature ladder, specified in one of several ways: a number of
    temperatures, a maximum temperature, or the inverse temperatures
    themselves. More temperatures move samples between modes more readily,
    at proportional cost.

``adaptive``, ``adaptation-lag``, ``adaptation-time``
    Let the temperature spacing adapt as the run proceeds.

``scale-factor``
    The stretch-move scale factor.

``logl-function``
    Which model attribute to use as the log likelihood.

epsie
=====

See the `epsie documentation <https://github.com/cdcapano/epsie>`_.

``nchains``
    Required. The number of parallel chains.

``ntemps``, ``inverse-temperatures-file``
    The temperature ladder: either a number of temperatures, or a file
    holding the inverse temperatures to use.

``swap-interval``
    Iterations between attempts to swap samples between temperatures.

``seed``
    Seed for epsie's own random state.

``logl-function``
    Which model attribute to use as the log likelihood.

Jump proposals are configured in their own
``[jump_proposal-{parameter}]`` sections, which is what distinguishes epsie
from the other MCMC samplers: the proposal for each parameter can be chosen
and tuned individually.

*******************
A worked comparison
*******************

:ref:`inference_example_samplers` runs the same simple analytic likelihood
through each sampler with a minimal configuration, and plots the resulting
posteriors on top of one another. It is a good starting point both for
copying a configuration and for checking that a sampler is installed and
behaving.
