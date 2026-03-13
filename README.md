# Persephone ABM

Agent-based model for simulating skeletal lesion formation and mortality in archaeological populations, based on the Usher3 illness-death model.

## Dependencies

- R (>= 4.0)
- demohaz
- testthat

## Running Tests

Source the module and run tests:

```bash
Rscript -e "source('R/mortality.R'); testthat::test_file('tests/testthat/unit/test-mortality.R')"
```

## Building Documentation

*Note: Requires R package structure with DESCRIPTION file.*

```bash
Rscript -e "roxygen2::roxygenize('.')"
```

