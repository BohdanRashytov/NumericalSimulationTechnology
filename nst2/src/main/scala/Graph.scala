package nst2

import java.awt.Color
import javax.swing.JFrame
import org.math.plot.Plot2DPanel

object Graph {
  val colors: Map[Int, Color] = Map(
    1 -> Color.red,
    2 -> Color.red.darker(),

    3 -> Color.pink,
    4 -> Color.pink.darker(),

    5 -> Color.orange,
    6 -> Color.orange.darker(),

    7 -> Color.yellow,
    8 -> Color.yellow.darker(),

    9 -> Color.green,
    10 -> Color.green.darker(),

    11 -> Color.magenta,
    12 -> Color.magenta.darker(),

    13 -> Color.cyan,
    14 -> Color.cyan.darker()
  )

  val countOfColors = colors.keySet.size

  def paintLines(f: (Double, Double) => (Double, Double), ps: List[(Double, Double)], defaultPoint: List[(Double, Double)], vihPoint: List[(Double, Double)], plotName: String) = {
    val points = ps.filterNot(pp => pp._1 == 2.0 && pp._2 >= 1.0 && pp._2 <= 2.0)
    val plot: Plot2DPanel = new Plot2DPanel()
    val frame: JFrame = new JFrame(plotName)
    frame.setSize(900, 900)
    frame.setVisible(true)
    frame.setContentPane(plot)

    val lines = points.map(point => Array(Array(point._1, point._2), {
      val xx = point._1 + f(point._1, point._2)._1*r(point._1, point._2)*0.01
      val yy = point._2 + f(point._1, point._2)._2*r(point._1, point._2)*0.01
      Array(xx, yy)
    }))
    lines.foreach(arr => plot.addLinePlot("lines", Color.BLACK, arr))
    plot.addScatterPlot(s"points$plotName", Color.BLUE, defaultPoint.map(p => Array(p._1, p._2)).toArray)
    plot.addScatterPlot(s"points$plotName", Color.GREEN, vihPoint.drop(0).map(p => Array(p._1, p._2)).toArray)
  }

  def r(x: Double, y: Double) = Math.sqrt(x*x + y*y)

  def paintPlot(f: (Double, Double) => Double, points: List[(Double, Double)], defaultPoint: List[(Double, Double)], vihPoint: List[(Double, Double)], plotName: String, ff: Option[(Double, Double) => (Double, Double)]) = {

    val plot: Plot2DPanel = new Plot2DPanel()
    val frame: JFrame = new JFrame(plotName)
    frame.setSize(900, 900)
    frame.setVisible(true)
    frame.setContentPane(plot)

    val max = {
      val np = points.map(p => f(p._1, p._2))
      var max = np.head
      np.foreach(q => if (q >= max) max = q)
      max
    }
    val min = {
      val np = points.map(p => f(p._1, p._2))
      var min = np.head
      np.foreach(q => if (q <= min) min = q)
      min
    }


    val groups: Map[Int, (List[(Double, Double)])] = {
      var initGroups: Map[Int, (List[(Double, Double)])] = Map().withDefaultValue(List[(Double, Double)]())
      val delta = (max - min + 0.00001) / countOfColors
      points.foreach(point => {
        val i = ((f(point._1, point._2) - min) / delta).toInt + 1
        initGroups += i -> (initGroups(i) ::: List((point._1, point._2)))
      })
      initGroups
    }
    groups.foreach(group => {
      val arr = group._2.map(point => Array(point._1, point._2)).toArray
      if (arr.size != 0) plot.addScatterPlot(s"$group", colors(group._1), arr)
    })

    ff.foreach(f => {
      val lines = points.map(point => Array(Array(point._1, point._2), {
        val xx = point._1 + f(point._1, point._2)._1*0.025
        val yy = point._2 + f(point._1, point._2)._2*0.025
        Array(xx, yy)
      }))
      lines.foreach(arr => plot.addLinePlot("lines", Color.BLACK, arr))
    })

    plot.addScatterPlot(s"points$plotName", Color.BLUE, defaultPoint.map(p => Array(p._1, p._2)).toArray)
    plot.addScatterPlot(s"points$plotName", Color.GREEN, vihPoint.drop(0).map(p => Array(p._1, p._2)).toArray)

    if (ff.isEmpty) {
      val scale: Plot2DPanel = new Plot2DPanel()
      val frameScale: JFrame = new JFrame("Scale for " + plotName)
      frameScale.setSize(900, 300)
      frameScale.setVisible(true)
      frameScale.setContentPane(scale)

      val K = 100
      val L = 100
      val delta = (max - min) / countOfColors

      (1 to countOfColors).foreach(i => {
        val ymin = min + (i - 1)*delta
        val ymax = min + i*delta
        val colorPoint = (for {
          k <- 1 until K
          l <- 1 until L
        } yield {
          Array(1.0*l*(ymax - ymin)/L + ymin, 1.0*k/K)
        }).toArray

        scale.addScatterPlot("scale", colors(i), colorPoint)
      })
    }
  }
}
