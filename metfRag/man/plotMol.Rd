\name{plotMol}
\alias{plotMol}
\title{Plotting a molecule}
\usage{
  plotMol(mol = NULL, smiles = NULL, watermark = NULL)
}
\arguments{
  \item{mol}{An rCDK \code{mol} object}

  \item{smiles}{A SMILES string, ignored if mol is
  supplied}

  \item{watermark}{A string written across the image as
  watermark}
}
\description{
  rCDK provides either functions to plot into a Java
  window, or into a rasterImage object. plotMol is a
  high-level plotting function to plot into an R device.
}

