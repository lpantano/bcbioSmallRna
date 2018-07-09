# bcbioSmallRna 0.0.1

## Fix
* Add plots to the template and remove it from the function
* Fix bcbioSmallSizeDist to avoid duplication of columns when getting adapter information

## Feature
* Change to use varianzeStabilization from DESeq2 to normalize counts,
  as consequence, name of slot of normalized count changed to `log`
* Load final folder using project folder path instead
* Feature: Improve cluster annotation