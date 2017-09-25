package nst1

import java.lang.Math._
import Jama.Matrix

object NST1 {

  val M: Int = 200
  val alpha: Double = 0.0
  val vInf: (Double, Double) = (cos(alpha), sin(alpha))
  val delta: Double = 0.001
  val gamma0: Double = 0.0

  val points: Map[Int, (Double, Double)] = {
    val n1 = M / 3
    val n2 = M - n1 + 1
    val list: List[(Double, Double)] =
      (1 to n1).map(i => (1.5 + (i - 1) * 0.5 / (n1 - 1), 1.5 + (i - 1) * 0.5 / (n1 - 1))).toList :::
        (2 to n2).map(i => (2.0, 2.0 - 1.0 * (i - 1) / (n2 - 1))).toList
    list.zipWithIndex.map(e => e._2 + 1 -> e._1).toMap
  }

  val colPoints: Map[Int, (Double, Double)] = (1 until M).map(i => i -> ((points(i)._1 + points(i + 1)._1) / 2, (points(i)._2 + points(i + 1)._2) / 2)).toMap

  val Rj: (Double, Double, Int) => Double = (x: Double, y: Double, j: Int) => max(delta, sqrt((x - points(j)._1) * (x - points(j)._1) + (y - points(j)._2) * (y - points(j)._2)))
  val Vj: (Double, Double, Int) => (Double, Double) = (x: Double, y: Double, j: Int) => ((points(j)._2 - y) / (2 * PI * Rj(x, y, j) * Rj(x, y, j)), (x - points(j)._1) / (2 * PI * Rj(x, y, j) * Rj(x, y, j)))


  val normal: Map[Int, (Double, Double)] = (1 until M).map(k => k -> {
    val deltaK = sqrt((points(k+1)._1 - points(k)._1)*(points(k+1)._1 - points(k)._1) + (points(k+1)._2 - points(k)._2)*(points(k+1)._2 - points(k)._2))
    (-(points(k+1)._2 - points(k)._2)/deltaK, (points(k+1)._1 - points(k)._1)/deltaK)
  }).toMap

  val gamma: Map[Int, Double] = {
    val b: List[Double] = colPoints.toList.sortBy(_._1).map(ixy => -mul(vInf, normal(ixy._1))) ::: List(gamma0)
    val a: Map[(Int, Int), Double] = (for {
      k <- 1 until M
      j <- 1 to M
    } yield (k, j) -> mul(Vj(colPoints(k)._1, colPoints(k)._2, j), normal(k))).toMap ++
      (for {
        j <- 1 to M
      } yield (M, j) -> 1.0).toMap
    val B = new Matrix(b.toArray, b.size)
    val A = {
      val array: Array[Array[Double]] = (1 to M).map(k => {
        (1 to M). map(j => a(k,j)).toArray
      }).toArray
      new Matrix(array)
    }
    A.solve(B).getRowPackedCopy.zipWithIndex.map(e => e._2 + 1 -> e._1).toMap
  }

  val v: (Double, Double) => (Double, Double) = (x: Double, y: Double) => (vInf._1 + (1 to M).map(i => gamma(i)*Vj(x, y, i)._1).sum, vInf._2 + (1 to M).map(i => gamma(i)*Vj(x, y, i)._2).sum)

  val phi: (Double, Double) => Double = (x: Double, y: Double) => x*vInf._1 + y*vInf._2 + (1 to M).map(i => gamma(i)*atan((y - points(i)._2)/(x - points(i)._1))).sum /(2*PI)
  val psi: (Double, Double) => Double = (x: Double, y: Double) => y*vInf._1 - x*vInf._2 - (1 to M).map(i => gamma(i)*log(Rj(x, y, i))).sum /(2*PI)

  val cp: (Double, Double) => Double = (x: Double, y: Double) => 1 - mul(v(x, y), v(x, y))/mul(vInf, vInf)

  def mul(a: (Double, Double), b: (Double, Double)) = a._1 * b._1 + a._2 * b._2
  def mul(a: Double, b: (Double, Double)) = (a*b._1, a*b._2)
  def sum(a: (Double, Double), b: (Double, Double)) = (a._1 + b._1, a._2 + b._2)



  val (x1, x2, y1, y2) = (0.05, 2.95, 0.05, 2.95)
  val h: Double = 0.05
  val n: Int = ((x2 - x1)/h).toInt

  val grid: List[Array[Array[Double]]] = (0 to n).map(i => List(Array(Array(i*h + x1, y1), Array(i*h + x1, y2)), Array(Array(x1, i*h + y1), Array(x2, i*h + y1)))).toList.flatten

  val allPoints: List[(Double, Double)] = (for {
    x <- 0 to n
    y <- 0 to n
  } yield (x1 + x*h, y1 + y*h)).toList


  def main(args: Array[String]): Unit = {
    Graph.paintPlot(psi, allPoints, points.values.toList, "PSI", None)
    Graph.paintPlot(phi, allPoints, points.values.toList, "PHI", None)
    Graph.paintPlot(cp, allPoints, points.values.toList, "CP", None)

    Graph.paintPlot(psi, allPoints, points.values.toList, "PSI - V", Some(v))
    Graph.paintPlot(phi, allPoints, points.values.toList, "PHI - V", Some(v))
    Graph.paintPlot(cp, allPoints, points.values.toList, "CP - V", Some(v))

    Graph.paintLines(v, allPoints, points.values.toList, "V")
  }
}
