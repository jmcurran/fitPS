# Stage 9.1.1 version-number repair

Stage 9.1 intended to begin the Stage 9 release series at `1.1.1.001`, but its generated runner formatted the four-part version incorrectly.

The faulty expression combined a three-element release-series vector with the formatted build number using both `sep` and `collapse`:

```r
paste(newValues[1:3], sprintf("%03d", newValues[4L]), sep = ".", collapse = ".")
```

R recycled the build component across the three release components, producing `1.001.1.001.1.001`.

Stage 9.1.1 repairs the committed metadata and uses `1.1.1.002`, because the Stage 9.1 attempt consumed build `001` under the fitPS attempted-build version policy.

Future stage runners must construct four-part versions by formatting each component first and collapsing a single four-element vector, for example:

```r
newVersion = paste(
  c(as.character(newValues[1:3]), sprintf("%03d", newValues[4L])),
  collapse = "."
)
```

No installed R implementation changes are made by this repair.
