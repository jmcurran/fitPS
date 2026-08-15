# Stage 10.4.2 build-path capture repair

Stage 10.4.1 completed source validation, committed the Stage 10.4 API rename, created its completed repository archive, and successfully built the source package. The runner then failed because shell command substitution captured the verbose `devtools::build()` console transcript together with the returned archive path.

Stage 10.4.2 preserves the committed Stage 10.4.1 source state and repairs only the stage completion workflow. The runner records the path returned by `devtools::build()` directly from R into a stage-specific path file, verifies that exact archive, installs it, and then performs portable cleanup of older completed archives.
