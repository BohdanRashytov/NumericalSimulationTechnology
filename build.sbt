name := "NumericalSimulationTechnology"

version := "1.0"

scalaVersion := "2.12.2"
libraryDependencies += "org.scala-lang" % "scala-library" % "2.12.2"

lazy val root = (project in file("."))
  .settings(name := "root")

lazy val nst1 = (project in file ("nst1"))
.settings(
  libraryDependencies += "gov.nist.math" % "jama" % "1.0.3",
  libraryDependencies += "com.github.yannrichet" % "JMathPlot" % "1.0",
  libraryDependencies += "com.github.yannrichet" % "JMathIO" % "1.0",
  libraryDependencies += "com.github.yannrichet" % "JMathArray" % "1.0"
)

lazy val nst2 = (project in file ("nst2"))
  .settings(
    libraryDependencies += "gov.nist.math" % "jama" % "1.0.3",
    libraryDependencies += "com.github.yannrichet" % "JMathPlot" % "1.0",
    libraryDependencies += "com.github.yannrichet" % "JMathIO" % "1.0",
    libraryDependencies += "com.github.yannrichet" % "JMathArray" % "1.0"
  )

lazy val nst3 = (project in file ("nst3"))
  .settings(
    libraryDependencies += "gov.nist.math" % "jama" % "1.0.3",
    libraryDependencies += "com.github.yannrichet" % "JMathPlot" % "1.0",
    libraryDependencies += "com.github.yannrichet" % "JMathIO" % "1.0",
    libraryDependencies += "com.github.yannrichet" % "JMathArray" % "1.0"
  )