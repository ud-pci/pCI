Working at UD
=============

The University of Delaware currently houses and maintains the Caviness community cluster and the NSF funded DARWIN cluster. 

.. list-table:: 
   :widths: 25 50
   :header-rows: 1

   * - Cluster
     - Processor
   * - Caviness
     - 2x18-core 2.10 GHz Intel E5-2695 v4 "Broadwell"
   * - DARWIN
     - 2x32-core AMD EPYC 7002 Series Processors

Here's more information about the `UD clusters <https://docs.hpc.udel.edu/>`_.

Compiling at UD
---------------

The computers at UD utilize the VALET system for installing software. We load the cmake and openmpi modules into your environment using the command ``vpkg_require``. 

.. code-block::

   $ vpkg_require intel
   Adding package `intel/2018u4` to your environment


To see what packages have been added to your environment, you can use the ``vpkg_history`` command.

.. code-block:: 

   $ vpkg_history
   [standard]
     intel/2018u4

To remove the changes produced by ``vpkg_require``, you can use the ``vpkg_rollback`` command.

.. code-block:: 

   $ vpkg_rollback


The parallel codes are built using the ``CMakeLists.txt`` file. 

A Debug build can be done:

.. code-block::

   $ ls
   pCI
   $ cd pCI
   $ mkdir build-debug
   $ cd build-debug
   $ vpkg_require cmake openmpi/4.1.0:intel-2020
   $ FC=mpifort cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=$(pwd)/../20200402-debug ..
      :
   $ make
   $ make install

An optimized build demands a little more:

.. code-block:: 

   $ cd ..
   $ mkdir build-opt
   $ cd build-opt
   $ FC=mpifort cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(pwd)/../20200317-opt  -DCMAKE_Fortran_FLAGS_RELEASE="-g -O3 -mcmodel=large -xHost -m64" ..
      :
   $ make
   $ make install


.. note::
    
   Each community cluster has their own versions of cmake and openmpi. To see what versions and variants of a package are available on the cluster, use the ``vpkg_versions`` command.

Running Jobs at UD
------------------

The UD clusters utilize the Slurm workload manager (job scheduling system) to manage and control the resources available to computational tasks. Users can either run their jobs interactively or submit in batch. For an interacive job, the user must type the commands they wish to execute in real time. For a batch job, those sequence of commands are saved to a file, known as a job script, which is submitted to Slurm. Using batch job scripts have several advantages such as reusability and increased job throughput. 

Here's more information about `running jobs <http://docs.hpc.udel.edu/abstract/darwin/runjobs/runjobs>`_.  

Job scripts
~~~~~~~~~~~

It is strongly recommended to use a job script file patterned after the prototypes in ``/opt/templates/``. There are ``README.md`` files in each subdirectory to explain the use of the templates. The following job scripts are valid for usage only in the Caviness and DARWIN community clusters. 

Serial job script
#################
See this :download:`serial job script template <src/serial.qs>` with full descriptions.

   .. raw:: html

      <details>
      <summary> 
         Click here to see the serial job script template..
      </summary>

   .. code-block:: 
   
      #!/bin/bash -l
      #
      # Sections of this script that can/should be edited are delimited by a
      # [EDIT] tag.  All Slurm job options are denoted by a line that starts
      # with "#SBATCH " followed by flags that would otherwise be passed on
      # the command line.  Slurm job options can easily be disabled in a
      # script by inserting a space in the prefix, e.g. "# SLURM " and
      # reenabled by deleting that space.
      #
      # This is a batch job template for a program using a single processor
      # core/thread (a serial job).
      #
      #SBATCH --ntasks=1
      #
      # [EDIT] All jobs have memory limits imposed.  The default is 1 GB per
      #        CPU allocated to the job.  The default can be overridden either
      #        with a per-node value (--mem) or a per-CPU value (--mem-per-cpu)
      #        with unitless values in MB and the suffixes K|M|G|T denoting
      #        kibi, mebi, gibi, and tebibyte units.  Delete the space between
      #        the "#" and the word SBATCH to enable one of them:
      #
      # SBATCH --mem=8G
      # SBATCH --mem-per-cpu=1024M
      #
      # [EDIT] Each node in the cluster has local scratch disk of some sort
      #        that is always mounted as /tmp.  Per-job and per-step temporary
      #        directories are automatically created and destroyed by the
      #        auto_tmpdir plugin in the /tmp filesystem.  To ensure a minimum
      #        amount of free space on /tmp when your job is scheduled, the
      #        --tmp option can be used; it has the same behavior unit-wise as
      #        --mem and --mem-per-cpu.  Delete the space between the "#" and the
      #        word SBATCH to enable:
      #
      # SBATCH --tmp=24G
      #
      # [EDIT] It can be helpful to provide a descriptive (terse) name for
      #        the job (be sure to use quotes if there's whitespace in the
      #        name):
      #
      #SBATCH --job-name=serial_job
      #
      # [EDIT] The partition determines which nodes can be used and with what
      #        maximum runtime limits, etc.  Partition limits can be displayed
      #        with the "sinfo --summarize" command.
      #
      # SBATCH --partition=standard
      #
      #        To run with priority-access to resources owned by your workgroup,
      #        use the "_workgroup_" partition:
      #
      # SBATCH --partition=_workgroup_
      #
      # [EDIT] The maximum runtime for the job; a single integer is interpreted
      #        as a number of minutes, otherwise use the format
      #
      #          d-hh:mm:ss
      #
      #        Jobs default to the default runtime limit of the chosen partition
      #        if this option is omitted.
      #
      #SBATCH --time=0-02:00:00
      #
      #        You can also provide a minimum acceptable runtime so the scheduler
      #        may be able to run your job sooner.  If you do not provide a
      #        value, it will be set to match the maximum runtime limit (discussed
      #        above).
      #
      # SBATCH --time-min=0-01:00:00
      #
      # [EDIT] By default SLURM sends the job's stdout to the file "slurm-<jobid>.out"
      #        and the job's stderr to the file "slurm-<jobid>.err" in the working
      #        directory.  Override by deleting the space between the "#" and the
      #        word SBATCH on the following lines; see the man page for sbatch for
      #        special tokens that can be used in the filenames:
      #
      # SBATCH --output=%x-%j.out
      # SBATCH --error=%x-%j.out
      #
      # [EDIT] Slurm can send emails to you when a job transitions through various
      #        states: NONE, BEGIN, END, FAIL, REQUEUE, ALL, TIME_LIMIT,
      #        TIME_LIMIT_50, TIME_LIMIT_80, TIME_LIMIT_90, ARRAY_TASKS.  One or more
      #        of these flags (separated by commas) are permissible for the
      #        --mail-type flag.  You MUST set your mail address using --mail-user
      #        for messages to get off the cluster.
      #
      # SBATCH --mail-user='my_address@udel.edu'
      # SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
      #
      # [EDIT] By default we DO NOT want to send the job submission environment
      #        to the compute node when the job runs.
      #
      #SBATCH --export=NONE
      #
   
      #
      # [EDIT] Define a Bash function and set this variable to its
      #        name if you want to have the function called when the
      #        job terminates (time limit reached or job preempted).
      #
      #        PLEASE NOTE:  when using a signal-handling Bash
      #        function, any long-running commands should be prefixed
      #        with UD_EXEC, e.g.
      #
      #                 UD_EXEC mpirun vasp
      #
      #        If you do not use UD_EXEC, then the signals will not
      #        get handled by the job shell!
      #
      #job_exit_handler() {
      #  # Copy all our output files back to the original job directory:
      #  cp * "$SLURM_SUBMIT_DIR"
      #
      #  # Don't call again on EXIT signal, please:
      #  trap - EXIT
      #  exit 0
      #}
      #export UD_JOB_EXIT_FN=job_exit_handler
   
      #
      # [EDIT] By default, the function defined above is registered
      #        to respond to the SIGTERM signal that Slurm sends
      #        when jobs reach their runtime limit or are
      #        preempted.  You can override with your own list of
      #        signals using this variable -- as in this example,
      #        which registers for both SIGTERM and the EXIT
      #        pseudo-signal that Bash sends when the script ends.
      #        In effect, no matter whether the job is terminated
      #        or completes, the UD_JOB_EXIT_FN will be called.
      #
      #export UD_JOB_EXIT_FN_SIGNALS="SIGTERM EXIT"
   
      #
      # If you have VALET packages to load into the job environment,
      # uncomment and edit the following line:
      #
      #vpkg_require intel/2019
   
      #
      # Do general job environment setup:
      #
      . /opt/shared/slurm/templates/libexec/common.sh
   
      #
      # [EDIT] Add your script statements hereafter, or execute a script or program
      #        using the srun command.
      #
      srun date <code>

.. raw:: html

   </details>

Once the job script has been set up, you can submit the job using the ``sbatch`` command:

.. code-block:: 

   sbatch serial.qs

Parallel job script
###################
See this :download:`parallel job script template <src/openmpi.qs>` with full descriptions.

   .. raw:: html

      <details>
      <summary>
         Click here to see the parallel job script template..
      </summary>

   .. code-block:: 

      #!/bin/bash -l
      #
      # Sections of this script that can/should be edited are delimited by a
      # [EDIT] tag.  All Slurm job options are denoted by a line that starts
      # with "#SBATCH " followed by flags that would otherwise be passed on
      # the command line.  Slurm job options can easily be disabled in a
      # script by inserting a space in the prefix, e.g. "# SLURM " and
      # reenabled by deleting that space.
      #
      # This is a batch job template for a program using multiple processor
      # cores/threads on one or more nodes.  This particular variant should
      # be used with Open MPI or another MPI library that is tightly-
      # integrated with Slurm.
      #
      # [EDIT] There are several ways to communicate the number and layout
      #        of worker processes.  Under GridEngine, the only option was
      #        to request a number of slots and GridEngine would spread the
      #        slots across an arbitrary number of nodes (not necessarily
      #        with a common number of worker per node, either).  This method
      #        is still permissible under Slurm by providing ONLY the
      #        --ntasks option:
      #
      #             #SBATCH --ntasks=<nproc>
      #
      #        To limit the number of nodes used to satisfy the distribution
      #        of <nproc> workers, the --nodes option can be used in addition
      #        to --ntasks:
      #
      #             #SBATCH --nodes=<nhosts>
      #             #SBATCH --ntasks=<nproc>
      #
      #        in which case, <nproc> workers will be allocated to <nhosts>
      #        nodes in round-robin fashion.
      #
      #        For a uniform distribution of workers the --tasks-per-node
      #        option should be used with the --nodes option:
      #
      #             #SBATCH --nodes=<nhosts>
      #             #SBATCH --tasks-per-node=<nproc-per-node>
      #
      #        The --ntasks option can be omitted in this case and will be
      #        implicitly equal to <nhosts> * <nproc-per-node>.
      #
      #        Given the above information, set the options you want to use
      #        and add a space between the "#" and the word SBATCH for the ones
      #        you don't want to use.
      #
      #SBATCH --nodes=<nhosts>
      #SBATCH --ntasks=<nproc>
      #SBATCH --tasks-per-node=<nproc-per-node>
      #
      # [EDIT] Normally, each MPI worker will not be multithreaded; if each
      #        worker allows thread parallelism, then alter this value to
      #        reflect how many threads each worker process will spawn.
      #
      #SBATCH --cpus-per-task=1
      #
      # [EDIT] All jobs have memory limits imposed.  The default is 1 GB per
      #        CPU allocated to the job.  The default can be overridden either
      #        with a per-node value (--mem) or a per-CPU value (--mem-per-cpu)
      #        with unitless values in MB and the suffixes K|M|G|T denoting
      #        kibi, mebi, gibi, and tebibyte units.  Delete the space between
      #        the "#" and the word SBATCH to enable one of them:
      #
      # SBATCH --mem=8G
      # SBATCH --mem-per-cpu=1024M
      #
      # [EDIT] Each node in the cluster has local scratch disk of some sort
      #        that is always mounted as /tmp.  Per-job and per-step temporary
      #        directories are automatically created and destroyed by the
      #        auto_tmpdir plugin in the /tmp filesystem.  To ensure a minimum
      #        amount of free space on /tmp when your job is scheduled, the
      #        --tmp option can be used; it has the same behavior unit-wise as
      #        --mem and --mem-per-cpu.  Delete the space between the "#" and the
      #        word SBATCH to enable:
      #
      # SBATCH --tmp=24G
      #
      # [EDIT] It can be helpful to provide a descriptive (terse) name for
      #        the job (be sure to use quotes if there's whitespace in the
      #        name):
      #
      #SBATCH --job-name=openmpi_job
      #
      # [EDIT] The partition determines which nodes can be used and with what
      #        maximum runtime limits, etc.  Partition limits can be displayed
      #        with the "sinfo --summarize" command.
      #
      # SBATCH --partition=standard
      #
      #        To run with priority-access to resources owned by your workgroup,
      #        use the "_workgroup_" partition:
      #
      # SBATCH --partition=_workgroup_
      #
      # [EDIT] The maximum runtime for the job; a single integer is interpreted
      #        as a number of minutes, otherwise use the format
      #
      #          d-hh:mm:ss
      #
      #        Jobs default to the default runtime limit of the chosen partition
      #        if this option is omitted.
      #
      #SBATCH --time=0-02:00:00
      #
      #        You can also provide a minimum acceptable runtime so the scheduler
      #        may be able to run your job sooner.  If you do not provide a
      #        value, it will be set to match the maximum runtime limit (discussed
      #        above).
      #
      # SBATCH --time-min=0-01:00:00
      #
      # [EDIT] By default SLURM sends the job's stdout to the file "slurm-<jobid>.out"
      #        and the job's stderr to the file "slurm-<jobid>.err" in the working
      #        directory.  Override by deleting the space between the "#" and the
      #        word SBATCH on the following lines; see the man page for sbatch for
      #        special tokens that can be used in the filenames:
      #
      # SBATCH --output=%x-%j.out
      # SBATCH --error=%x-%j.out
      #
      # [EDIT] Slurm can send emails to you when a job transitions through various
      #        states: NONE, BEGIN, END, FAIL, REQUEUE, ALL, TIME_LIMIT,
      #        TIME_LIMIT_50, TIME_LIMIT_80, TIME_LIMIT_90, ARRAY_TASKS.  One or more
      #        of these flags (separated by commas) are permissible for the
      #        --mail-type flag.  You MUST set your mail address using --mail-user
      #        for messages to get off the cluster.
      #
      # SBATCH --mail-user='my_address@udel.edu'
      # SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
      #
      # [EDIT] By default we DO NOT want to send the job submission environment
      #        to the compute node when the job runs.
      #
      #SBATCH --export=NONE
      #

      #
      # [EDIT] Do any pre-processing, staging, environment setup with VALET
      #        or explicit changes to PATH, LD_LIBRARY_PATH, etc.
      #
      vpkg_require openmpi/default

      #
      # [EDIT] If you're not interested in how the job environment gets setup,
      #        uncomment the following.
      #
      #UD_QUIET_JOB_SETUP=YES

      #
      # [EDIT] Slurm has a specific MPI-launch mechanism in srun that can speed-up
      #        the startup of jobs with large node/worker counts.  Uncomment this
      #        line if you want to use that in lieu of mpirun.
      #
      #UD_USE_SRUN_LAUNCHER=YES

      #
      # [EDIT] By default each MPI worker process will be bound to a core/thread
      #        for better efficiency.  Uncomment this to NOT affect such binding.
      #
      #UD_DISABLE_CPU_AFFINITY=YES

      #
      # [EDIT] MPI ranks are distributed <nodename>(<rank>:<socket>.<core>,..)
      #
      #    CORE    sequentially to all allocated cores on each allocated node in
      #            the sequence they occur in SLURM_NODELIST (this is the default)
      #
      #              -N2 -n4 => n000(0:0.0,1:0.1,2:0.2,3:0.3); n001(4:0.0,5:0.1,6:0.2,7:0.3)
      #
      #    NODE    round-robin across the nodes allocated to the job in the sequence
      #            they occur in SLURM_NODELIST
      #
      #              -N2 -n4 => n000(0:0.0,2:0.1,4:0.2,6:0.3); n001(1:0.0,3:0.1,5:0.2,7:0.3)
      #
      #    SOCKET  round-robin across the allocated sockets on each allocated node
      #            in the sequence they occur in SLURM_NODELIST
      #
      #              -N2 -n4 => n000(0:0.0,2:0.1,4:1.0,6:1.1); n001(1:0.0,3:0.1,5:1.0,7:1.1)
      #
      #            PLEASE NOTE:  socket mode requires use of the --exclusive flag
      #            to ensure uniform allocation of cores across sockets!
      #
      #UD_MPI_RANK_DISTRIB_BY=CORE

      #
      # [EDIT] By default all MPI byte transfers are limited to NOT use any
      #        TCP interfaces on the system.  Setting this variable will force
      #        the job to NOT use any Infiniband interfaces.
      #
      #UD_DISABLE_IB_INTERFACES=YES

      #
      # [EDIT] Should Open MPI display LOTS of debugging information as the job
      #        executes?  Uncomment to enable.
      #
      #UD_SHOW_MPI_DEBUGGING=YES

      #
      # [EDIT] Define a Bash function and set this variable to its
      #        name if you want to have the function called when the
      #        job terminates (time limit reached or job preempted).
      #
      #        PLEASE NOTE:  when using a signal-handling Bash
      #        function, any long-running commands should be prefixed
      #        with UD_EXEC, e.g.
      #
      #                 UD_EXEC mpirun vasp
      #
      #        If you do not use UD_EXEC, then the signals will not
      #        get handled by the job shell!
      #
      #job_exit_handler() {
      #  # Copy all our output files back to the original job directory:
      #  cp * "$SLURM_SUBMIT_DIR"
      #
      #  # Don't call again on EXIT signal, please:
      #  trap - EXIT
      #  exit 0
      #}
      #export UD_JOB_EXIT_FN=job_exit_handler

      #
      # [EDIT] By default, the function defined above is registered
      #        to respond to the SIGTERM signal that Slurm sends
      #        when jobs reach their runtime limit or are
      #        preempted.  You can override with your own list of
      #        signals using this variable -- as in this example,
      #        which registers for both SIGTERM and the EXIT
      #        pseudo-signal that Bash sends when the script ends.
      #        In effect, no matter whether the job is terminated
      #        or completes, the UD_JOB_EXIT_FN will be called.
      #
      #export UD_JOB_EXIT_FN_SIGNALS="SIGTERM EXIT"

      #
      # Do standard Open MPI environment setup (networks, etc.)
      #
      . /opt/shared/slurm/templates/libexec/openmpi.sh

      #
      # [EDIT] Execute your MPI program
      #
      ${UD_MPIRUN} ./my_mpi_program arg1 "arg2 has spaces" arg3
      mpi_rc=$?

      #
      # [EDIT] Do any cleanup work here...
      #

      #
      # Be sure to return the mpirun's result code:
      #
      exit $mpi_rc

.. raw:: html

   </details>


Once the job script has been set up, you can submit the job using the ``sbatch`` command:

.. code-block:: 

   sbatch openmpi.qs


Managing Jobs at UD
-------------------

Once the job has been submitted, you can monitor the status of your job using the ``squeue`` command:

.. code-block:: bash

   squeue -u <username>
   squeue -p <partition_name>

You can also continuously monitor your job by using the ``watch`` command:

.. code-block:: 

   watch squeue -u <username>
   watch squeue -p <partition_name>

(Caviness only) To see information about the current utilization of guaranteed resources for the workgroup, you can run the ``squota`` command:

.. code-block:: 

   squota

To cancel your job, you can run the ``scancel`` command:

.. code-block:: 

   scancel <job-id>

To see information about the partitions and nodes, you can run the ``sinfo`` command:

.. code-block:: 

   sinfo
   sinfo -p <partition-name>

To see information about your queued jobs, you can run the ``scontrol`` command:

.. code-block:: 

   scontrol show job <job-id>

To see information about a completed job, you can run the ``sacct`` command:

.. code-block:: 

   sacct -j <job-id>

More information about managing jobs can be found here for `Caviness <http://docs.hpc.udel.edu/abstract/caviness/runjobs/job_status>`_ and `DARWIN <http://docs.hpc.udel.edu/abstract/darwin/runjobs/job_status>`_.  

Additional Information and Support
----------------------------------
For additional information and support for running jobs on the UD clusters, please visit the respective cluster documentation pages:

| Caviness: `http://docs.hpc.udel.edu/abstract/caviness/caviness <http://docs.hpc.udel.edu/abstract/caviness/caviness>`_
| DARWIN: `http://docs.hpc.udel.edu/abstract/darwin/darwin <http://docs.hpc.udel.edu/abstract/darwin/darwin>`_
