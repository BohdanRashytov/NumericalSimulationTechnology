package nst2

import java.lang.Math._
import Jama.Matrix

object NST2 {

  val M: Int = 30
  val alpha: Double = 0.0
  val vInf: (Double, Double) = (cos(alpha), sin(alpha))
  val delta: Map[Int, Double] = Map[Int, Double]().withDefaultValue(0.001)
  val gamma0 = 0.0
  val dt = 3.0 / M
  val p = 3

  val points: Map[Int, (Double, Double)] = {
    val n1 = M / 3
    val n2 = M - n1 + 1
    val list: List[(Double, Double)] =
      (1 to n1).map(i => (1.5 + (i - 1) * 0.5 / (n1 - 1), 1.5 + (i - 1) * 0.5 / (n1 - 1))).toList :::
        (2 to n2).map(i => (2.0, 2.0 - 1.0 * (i - 1) / (n2 - 1))).toList
    list.zipWithIndex.map(e => e._2 + 1 -> e._1).toMap
  }

  val getPointIndex: Map[Int, Int] = Map(1 -> 1, 2 -> (M / 3), 3 -> M)
  var pointsP: Map[Int, Map[Int, (Double, Double)]] = Map(0 -> Map(1 -> points(getPointIndex(1)), 2 -> points(getPointIndex(2)), 3 -> points(getPointIndex(3))))

  val colPoints: Map[Int, (Double, Double)] = (1 until M).map(i => i -> ((points(i)._1 + points(i + 1)._1) / 2, (points(i)._2 + points(i + 1)._2) / 2)).toMap

  val Rj: (Double, Double, Int) => Double = (x: Double, y: Double, j: Int) => max(delta(j), sqrt((x - points(j)._1) * (x - points(j)._1) + (y - points(j)._2) * (y - points(j)._2)))

  val Rjp: (Double, Double, Int, Int) => Double = (x: Double, y: Double, j: Int, p: Int) => max(delta(j), sqrt((x - pointsP(j)(p)._1) * (x - pointsP(j)(p)._1) + (y - pointsP(j)(p)._2) * (y - pointsP(j)(p)._2)))

  val Vj: (Double, Double, Int) => (Double, Double) = (x: Double, y: Double, j: Int) =>
    ((points(j)._2 - y) / (2 * PI * Rj(x, y, j) * Rj(x, y, j)),
      (x - points(j)._1) / (2 * PI * Rj(x, y, j) * Rj(x, y, j)))

  val Vjp: (Double, Double, Int, Int) => (Double, Double) = (x: Double, y: Double, j: Int, p: Int) =>
    ((pointsP(j)(p)._2 - y) / (2 * PI * Rjp(x, y, j, p) * Rjp(x, y, j, p)),
      (x - pointsP(j)(p)._1) / (2 * PI * Rjp(x, y, j, p) * Rjp(x, y, j, p)))


  val normal: Map[Int, (Double, Double)] = (1 until M).map(k => k -> {
    val deltaK = sqrt((points(k + 1)._1 - points(k)._1) * (points(k + 1)._1 - points(k)._1) + (points(k + 1)._2 - points(k)._2) * (points(k + 1)._2 - points(k)._2))
    (-(points(k + 1)._2 - points(k)._2) / deltaK, (points(k + 1)._1 - points(k)._1) / deltaK)
  }).toMap

  def gammaCalc: Int => Map[Int, Double] = (t: Int) => {
    val b: List[Double] = colPoints.toList.sortBy(_._1).map(ixy => -mul(vInf, normal(ixy._1)) - (for {
      p <- 1 to p
      t <- 1 to t
    } yield gamma(t - 1)(getPointIndex(p)) * mul(Vjp(ixy._2._1, ixy._2._2, t, p), normal(ixy._1))).sum
    ) ::: List(gamma0 - (1 to p).map(p => (1 to t).map(t => gamma(t - 1)(getPointIndex(p))).sum).sum)
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
        (1 to M).map(j => a(k, j)).toArray
      }).toArray
      new Matrix(array)
    }
    A.solve(B).getRowPackedCopy.zipWithIndex.map(e => e._2 + 1 -> e._1).toMap
  }

  val v: Int => (Double, Double) => (Double, Double) = (t: Int) => (x: Double, y: Double) =>
    (vInf._1 + (1 to M).map(i => gamma(t)(i) * Vj(x, y, i)._1).sum + (1 to p).map(p => (1 to t).map(t => mul(gamma(t - 1)(getPointIndex(p)), Vjp(x, y, t, p))._1).sum).sum,
      vInf._2 + (1 to M).map(i => gamma(t)(i) * Vj(x, y, i)._2).sum + (1 to p).map(p => (1 to t).map(t => mul(gamma(t - 1)(getPointIndex(p)), Vjp(x, y, t, p))._2).sum).sum)

  var gamma: Map[Int, Map[Int, Double]] = Map[Int, Map[Int, Double]]()

  def findVP(point: (Double, Double), t: Int, index: Int): (Double, Double) = {
    val ddt = dt /  Graph.r(v(t)(point._1, point._2)._1, v(t)(point._1, point._2)._2)
    val xx = point._1 + v(t)(point._1, point._2)._1*ddt
    val yy = point._2 + v(t)(point._1, point._2)._2*ddt
    if (index == 1 && yy >= xx) return ((yy - xx) + yy, yy)
    if (index != 1 && xx <= 2.0) return ((2.0 - xx) + 2.0 , yy)
    (xx, yy)
  }

  (0 to 75).foreach(i => {
    if (i % 10 == 0) println(i)
    gamma += i -> gammaCalc(i)
    pointsP += i + 1 -> {
      Map(
        1 -> findVP((pointsP(i)(1)._1, pointsP(i)(1)._2), i, 1),
        2 -> findVP((pointsP(i)(2)._1, pointsP(i)(2)._2), i, 2),
        3 -> findVP((pointsP(i)(3)._1, pointsP(i)(3)._2), i, 3)
      )
    }
  })


  val phi: Int => (Double, Double) => Double = (t: Int) => (x: Double, y: Double) => x * vInf._1 + y * vInf._2 +
    (1 to M).map(i => gamma(t)(i) * atan((y - points(i)._2) / (x - points(i)._1))).sum / (2 * PI) +
    (1 to p).map(p => (1 to t).map(t => gamma(t - 1)(getPointIndex(p)) * atan((y - pointsP(t)(p)._2) / (x - pointsP(t)(p)._1))).sum).sum / (2 * PI)

  val psi: Int => (Double, Double) => Double = (t: Int) => (x: Double, y: Double) => y * vInf._1 - x * vInf._2 -
    (1 to M).map(i => gamma(t)(i) * log(Rj(x, y, i))).sum / (2 * PI) +
    (1 to p).map(p => (1 to t).map(t => gamma(t - 1)(getPointIndex(p)) * log(Rjp(x, y, t, p))).sum).sum / (2 * PI)

  def mul(a: (Double, Double), b: (Double, Double)) = a._1 * b._1 + a._2 * b._2

  def mul(a: Double, b: (Double, Double)) = (a * b._1, a * b._2)

  def sum(a: (Double, Double), b: (Double, Double)) = (a._1 + b._1, a._2 + b._2)


  val (x1, x2, y1, y2) = (0.05, 2.95, 0.05, 2.95)
  val h: Double = 0.05
  val n: Int = ((x2 - x1) / h).toInt

  val grid: List[Array[Array[Double]]] = (0 to n).map(i => List(Array(Array(i * h + x1, y1), Array(i * h + x1, y2)), Array(Array(x1, i * h + y1), Array(x2, i * h + y1)))).toList.flatten

  val allPoints: List[(Double, Double)] = (for {
    x <- 0 to n
    y <- 0 to n
  } yield (x1 + x * h, y1 + y * h)).toList


  def main(args: Array[String]): Unit = {
    val j = 25
    val i = 45
    val k = 65

    Graph.paintLines(v(j), allPoints, points.values.toList, pointsP.filter(_._1 <= j).values.flatten.toList.map(_._2), s"V-$j")
    Graph.paintLines(v(i), allPoints, points.values.toList, pointsP.filter(_._1 <= i).values.flatten.toList.map(_._2), s"V-$i")
    Graph.paintLines(v(k), allPoints, points.values.toList, pointsP.filter(_._1 <= k).values.flatten.toList.map(_._2), s"V-$k")

    Graph.paintPlot(psi(j), allPoints, points.values.toList, pointsP.filter(_._1 <= j).values.flatten.toList.map(_._2), s"PSI-$j", None)
    Graph.paintPlot(psi(i), allPoints, points.values.toList, pointsP.filter(_._1 <= i).values.flatten.toList.map(_._2), s"PSI-$i", None)
    Graph.paintPlot(psi(k), allPoints, points.values.toList, pointsP.filter(_._1 <= k).values.flatten.toList.map(_._2), s"PSI-$k", None)

    Graph.paintPlot(phi(j), allPoints, points.values.toList, pointsP.filter(_._1 <= j).values.flatten.toList.map(_._2), s"PHI-$j", None)
    Graph.paintPlot(phi(i), allPoints, points.values.toList, pointsP.filter(_._1 <= i).values.flatten.toList.map(_._2), s"PHI-$i", None)
    Graph.paintPlot(phi(k), allPoints, points.values.toList, pointsP.filter(_._1 <= k).values.flatten.toList.map(_._2), s"PHI-$k", None)

    Graph.paintPlot(psi(j), allPoints, points.values.toList, pointsP.filter(_._1 <= j).values.flatten.toList.map(_._2), s"PSI-V-$j", Some(v(j)))
    Graph.paintPlot(psi(i), allPoints, points.values.toList, pointsP.filter(_._1 <= i).values.flatten.toList.map(_._2), s"PSI-V-$i", Some(v(i)))
    Graph.paintPlot(psi(k), allPoints, points.values.toList, pointsP.filter(_._1 <= k).values.flatten.toList.map(_._2), s"PSI-V-$k", Some(v(k)))

    Graph.paintPlot(phi(j), allPoints, points.values.toList, pointsP.filter(_._1 <= j).values.flatten.toList.map(_._2), s"PHI-V-$j", Some(v(j)))
    Graph.paintPlot(phi(i), allPoints, points.values.toList, pointsP.filter(_._1 <= i).values.flatten.toList.map(_._2), s"PHI-V-$i", Some(v(i)))
    Graph.paintPlot(phi(k), allPoints, points.values.toList, pointsP.filter(_._1 <= k).values.flatten.toList.map(_._2), s"PHI-V-$k", Some(v(k)))
  }
}
