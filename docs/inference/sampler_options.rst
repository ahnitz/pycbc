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

The samplers fall broadly into two families:

**Nested samplers** integrate the likelihood over the prior. They return an
estimate of the evidence (:math:`\log Z`) as well as posterior samples, they
decide for themselves when they have converged, and they need no burn-in or
autocorrelation analysis. Most of the gravitational-wave configurations
under ``examples/inference`` use one.

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
     - What most of the ``examples/inference`` configurations use.
   * - ``ultranest``
     - nested
     - Has a step-sampling mode for higher dimension.
   * - ``nessai``
     - nested
     - Uses normalising flows to propose new points.
   * - ``multinest``
     - nested
     - Needs the external MultiNest library and ``pymultinest``.
   * - ``cpnest``
     - nested
     - Nested sampling parallelised over threads.
   * - ``snowline``
     - importance
     - Fits a Gaussian and importance samples from it.
   * - ``emcee``
     - MCMC
     - Affine-invariant ensemble sampler.
   * - ``ptemcee``
     - MCMC
     - Parallel-tempered ``emcee``.
   * - ``epsie``
     - MCMC
     - Parallel-tempered, with per-parameter jump proposals.

``pycbc.inference.sampler.samplers`` is the authoritative list; it also holds
a few special-purpose entries not covered here (``dummy``, which is selected
automatically when there are no variable parameters, plus ``refine`` and
``games``). ``emcee``, ``epsie``, ``ptemcee``, ``cpnest``, ``dynesty`` and
``nessai`` are registered only if their package imports, so with the package
missing the name is reported as unavailable rather than failing later.

**********************
Options for every run
**********************

``name``
    Required. Selects the sampler; see the table above.

Only ``emcee`` and ``dynesty`` are in ``requirements.txt``. The rest come
from ``companion.txt`` (``epsie``, ``cpnest``, ``pymultinest``,
``ultranest``, ``ptemcee``, ``nessai``, ``snowline``) and have to be
installed before the corresponding ``name`` will resolve.

Two more controls are worth knowing about even though they are *not* set in
the ``[sampler]`` section:

* ``--nprocesses`` on the ``pycbc_inference`` command line sets the number of
  parallel processes. It is passed through by ``dynesty``, ``nessai``,
  ``emcee``, ``ptemcee`` and ``epsie``. ``cpnest`` takes its own
  ``nthreads`` option instead, and ``ultranest``, ``snowline`` and
  ``multinest`` ignore it.
* ``--seed`` seeds numpy's global generator at startup (default 0), which is
  what sets the initial walker positions. Samplers that carry their own
  random state are not necessarily covered by it: ``epsie`` takes its own
  ``seed`` option, and ``dynesty`` keeps an internal ``rstate``.

*******************
Nested samplers
*******************

dynesty
=======

See the `dynesty documentation <https://dynesty.readthedocs.io/>`_, and the
:ref:`example configuration <platter_dynesty>`.

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

``maxiter``
    A hard cap on the number of iterations.

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
    bootstrapping. dynesty rejects the two together unless one of them is
    neutral, i.e. ``bootstrap = 0`` or ``enlarge = 1``.

Checkpointing
-------------

``checkpoint_time_interval``
    Seconds between checkpoints. If it is not set, the run does not
    checkpoint at all. It is also ignored when ``nlive`` is negative, since
    the dynamic nested sampler does not support checkpointing here.

``maxcall``
    How many likelihood calls dynesty runs for before returning so that
    pycbc can decide whether to checkpoint; it has no effect unless
    ``checkpoint_time_interval`` is set. Defaults to 5000 per pool process.

ultranest
=========

See the `UltraNest documentation
<https://johannesbuchner.github.io/UltraNest/>`_, and the
:ref:`example configuration <platter_ultranest>`.

``min_num_live_points``
    The number of live points.

The stopping conditions are ``dlogz``, ``dKL``, ``Lepsilon``, ``min_ess``,
``frac_remain``, ``max_iters``, ``max_ncalls`` and
``max_num_improvement_loops``.

``stepsampling``
    If true, use a hit-and-run step sampler instead of the default rejection
    sampling. UltraNest describes this as scaling better with
    dimensionality.

``update_interval_iter_fraction``, ``update_interval_ncall``
    How often the bounding region is updated.

``log_dir``, ``log_interval``, ``show_status``, ``enable_plots``
    Diagnostic output.

nessai
======

See the `nessai documentation <https://nessai.readthedocs.io/>`_, and the
:ref:`example configuration <platter_nessai>`.

``nlive``
    The number of live points.

``importance_nested_sampler``
    Selects nessai's importance nested sampler. PyCBC does not support it:
    ``from_config`` raises ``NotImplementedError`` for it, so leave this
    unset.

This section is validated against nessai itself rather than a fixed list:
any keyword accepted by ``FlowSampler`` or its ``run`` method may be given,
and anything else raises a ``RuntimeError`` naming the unknown options.

snowline
========

See the `snowline documentation
<https://johannesbuchner.github.io/snowline/>`_, and the
:ref:`example configuration <platter_snowline>`.

``num_global_samples``, ``num_gauss_samples``
    How many samples are drawn globally and from the fitted Gaussian.

``min_ess``
    Target effective sample size.

``max_ncalls``, ``max_improvement_loops``
    Stopping conditions.

multinest
=========

See the `PyMultiNest documentation
<https://johannesbuchner.github.io/PyMultiNest/>`_, and the
:ref:`example configuration <platter_multinest>`. This sampler needs the
MultiNest library itself as well as the ``pymultinest`` wrapper.

``nlivepoints``
    Required. The number of live points.

``evidence-tolerance``
    Target accuracy on the evidence; the main stopping condition.

``sampling-efficiency``
    Trades sampling efficiency against accuracy; see the MultiNest
    documentation for the recommended values.

``importance-nested-sampling``
    Turns on MultiNest's importance nested sampling.

``checkpoint-interval``
    Iterations between checkpoints.

Note that the joint prior cannot be used here, so constraints are read
separately from the ``[constraint-*]`` sections.

cpnest
======

See the `cpnest documentation <https://github.com/johnveitch/cpnest>`_, and
the :ref:`example configuration <platter_cpnest>`. All four options below
are required; there are no defaults.

``nlive``
    The number of live points.

``maxmcmc``
    The maximum length of the internal MCMC used to evolve a point.

``nthreads``
    Number of parallel threads cpnest should use. ``--nprocesses`` does not
    reach cpnest, so this is the knob to set.

``verbose``
    cpnest's verbosity level, as an integer.

*******************
MCMC samplers
*******************

Options shared by the MCMC samplers
===================================

``niterations``
    Run for a fixed number of iterations. Give either this or
    ``effective-nsamples``, not both.

``effective-nsamples``
    Alternatively, run until this many *effective* samples have been
    collected, i.e. after burn-in and thinning by the autocorrelation
    length. Use this rather than ``niterations`` if you want the run to
    decide for itself when it has enough. It requires
    ``checkpoint-interval`` to be set as well, since the target is tested at
    each checkpoint.

``checkpoint-interval``
    Iterations between checkpoints.

``checkpoint-signal``
    Signal to checkpoint and exit on, e.g. ``USR2``, for running under a
    batch system.

``thin-interval``
    Thin the stored chain by this factor as it is written. Must be smaller
    than ``checkpoint-interval``.

``max-samples-per-chain``
    Caps memory use by thinning on the fly once a chain reaches this length.
    Give either this or ``thin-interval``, not both.

Burn-in is configured in its own ``[sampler-burn_in]`` section:

``burn-in-test``
    Which test decides that the chain has burned in. The tests available to
    every MCMC sampler are ``halfchain``, ``min_iterations``,
    ``max_posterior``, ``posterior_step`` and ``nacl``; ``emcee`` adds
    ``ks_test``. They combine with ``&`` and ``|``, e.g.
    ``nacl & max_posterior``. If the section is absent, no burn-in test is
    applied. See :py:mod:`pycbc.inference.burn_in` for what each test means.

emcee
=====

See the `emcee documentation <https://emcee.readthedocs.io/>`_, and the
:ref:`example configuration <platter_emcee>`.

``nwalkers``
    Required. The number of walkers in the ensemble. emcee imposes a lower
    bound on this relative to the number of ``variable_params``; see its
    documentation.

``logpost-function``
    Which model attribute to use as the log posterior.

ptemcee
=======

See the `ptemcee documentation <https://github.com/willvousden/ptemcee>`_,
and the :ref:`example configuration <platter_ptemcee>`.

``nwalkers``
    Required, as for ``emcee``.

``ntemps``, ``tmax``, ``betas``, ``betas-file``
    The temperature ladder. Give a number of temperatures (``ntemps``), a
    maximum temperature (``tmax``), both, or the inverse temperatures
    themselves as a list (``betas``) or in an HDF file under
    ``.attrs['betas']`` (``betas-file``). ``betas`` and ``betas-file`` are
    mutually exclusive, and neither may be combined with ``ntemps`` or
    ``tmax``.

``adaptive``, ``adaptation-lag``, ``adaptation-time``
    Let the temperature spacing adapt as the run proceeds.

``scale-factor``
    The stretch-move scale factor.

``logl-function``
    Which model attribute to use as the log likelihood.

epsie
=====

See the `epsie documentation <https://github.com/cdcapano/epsie>`_, and the
:ref:`example configuration <platter_epsie>`.

``nchains``
    Required. The number of parallel chains.

``ntemps``, ``inverse-temperatures-file``
    The temperature ladder: either a number of temperatures, or an HDF file
    holding the inverse temperatures under ``.attrs['betas']``. Give one or
    the other, not both.

``swap-interval``
    Iterations between attempts to swap samples between temperatures.
    Default is 1.

``seed``
    Seed for epsie's own random number generator. If unset, epsie creates
    one, so ``--seed`` alone does not pin an epsie run.

``logl-function``
    Which model attribute to use as the log likelihood.

Jump proposals are configured in ``[jump_proposal-{params}]`` sections,
where ``{params}`` is a :py:const:`pycbc.VARARGS_DELIM` separated list. This
is what distinguishes epsie from the other MCMC samplers: the proposal for
each parameter can be chosen and tuned individually. A proposal must be
given for every *sampling* parameter (not every variable parameter); see
:py:func:`pycbc.inference.jump.epsie_proposals_from_config`.

*******************
A worked comparison
*******************

:ref:`inference_example_samplers` runs the same unimodal Gaussian likelihood
through most of these samplers with a minimal configuration and plots the
resulting posteriors on top of one another. It is a good starting point both
for copying a configuration and for checking that a sampler is installed and
behaving.
