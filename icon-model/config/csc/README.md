# Config wrappers in LUMI

The Buildbot builds and runs ICON inside a container.

An interactive shell session inside the container can be started with the
following command:

```bash
$ ./config/csc/exec.lumi.container.cce bash
```

where `exec.lumi.container.cce` is a symbolic link pointing to the current
default container wrapper.

Once the shell session is started, ICON can be configure and built as follows:

- for CPU:
    ```bash
    Apptainer> ./config/csc/lumi.cpu.cce
    ...
    Apptainer> make -j8
    ```
- for GPU:
    ```bash
    Apptainer> ./config/csc/lumi.gpu.cce
    ...
    Apptainer> make -j8
    ```

where `lumi.cpu.cce` and `lumi.gpu.cce` are symbolic links pointing the current
default configure wrappers (e.g. `lumi.cpu.container.cce-16.0.1.1` and
`lumi.gpu.container.cce-16.0.1.1`, respectively).

Alternatively, ICON can be configured, built and run as usual in the native
environment using the `lumi.*.native.*` configure wrappers.
