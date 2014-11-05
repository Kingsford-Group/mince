#!/bin/bash
exec scala -cp /home/robp/mince/scripts/parboiled-scala_2.11-1.1.6.jar:/home/robp/mince/scripts/spray-json_2.11-1.2.6.jar -nc -savecompiled "$0" "$@"
!#

import scala.sys.process._
import scala.collection.mutable.{ArrayBuffer, OpenHashMap => HMap}
import java.text.NumberFormat
import java.io.File

case class FileInfo(var name: String, var method: String, var size: Long, var isPE: Boolean = false ) {
}

def getFileList(f: File, method: String, extensions: List[String]): HMap[File, FileInfo] = {
    val pat = s"""(.*)(${extensions.mkString("|")})""".r
    println(s"""Pattern: ${pat}""")
    val retFiles = new HMap[File, FileInfo]
    f.listFiles.foreach( fn =>
            fn.toString() match {
                case pat(fname, ext) => {
                    retFiles.put(fn , FileInfo(fname, method, fn.length))
                }
                case _ => {}
            }
    )
    retFiles
}


val minceMap = new HMap[String, FileInfo]
val flipsMap = new HMap[String, FileInfo]
val minceNoRCMap = new HMap[String, FileInfo]
val scalceMap = new HMap[String, FileInfo]
val fastqzMap = new HMap[String, FileInfo]

val retMince = getFileList(new File(args(0)), "mince", List("offs.lz", "seqs.lz"))

retMince.foreach{ case (k,v) => 
    var base = k.getName.split("""\.""")(0)
    println(base)
    if (minceMap.contains(base)) {
        minceMap(base).size += v.size
    } else {
        minceMap.put(base, v)
    }
}

val retFlips = getFileList(new File(args(0)), "flips", List("flips.lz"))

retFlips.foreach{ case (k,v) =>
    var base = k.getName.split("""\.""")(0)
    if (flipsMap.contains(base)) {
        flipsMap(base).size += v.size
    } else {
        flipsMap.put(base, v)
    }   
}   

val retMinceNoRC = getFileList(new File(args(1)), "mince (no rc)", List("offs.lz", "seqs.lz"))

retMinceNoRC.foreach{ case (k,v) => 
    var base = k.getName.split("""\.""")(0)
    if (minceNoRCMap.contains(base)) {
        minceNoRCMap(base).size += v.size
    } else {
        minceNoRCMap.put(base, v)
    }
}


val retScalce = getFileList(new File(args(2)), "scalce", List("scalcer"))

retScalce.foreach{ case (k,v) => 
    var base = k.getName.split("_")(0)
    if (scalceMap.contains(base)) {
        scalceMap(base).size += v.size
        scalceMap(base).isPE = true
    } else {
        scalceMap.put(base, v)
    }
}

val retFastqz = getFileList(new File(args(3)), "fastqz", List("fxb.zpaq"))

retFastqz.foreach{ case (k,v) =>
    var base = k.getName.split("""\.""")(0)
    if (fastqzMap.contains(base)) {
        fastqzMap(base).size += v.size
    } else {
        fastqzMap.put(base, v)
    }
}

val mergedMap = (minceMap.toSeq ++ minceNoRCMap ++ scalceMap.toSeq ++ fastqzMap.toSeq).groupBy{ case (k,v) => k }

val numKeys = mergedMap.keySet.size
val methodsInOrder = Array("fastqz", "scalce", "mince", "mince (no rc)")

println("""\begin{center}""")
println("""\begin{tabular}{lcccc}""")

println("""\toprule""")
println(s"""Dataset & ${methodsInOrder.mkString(" & ")}\\\\""")
println("""\midrule""")

mergedMap.foreach{ case (k,v) =>
    val methodMap = v.map{ case (n,t) => 
        (t.method, if (t.size > 0) { s"""$$${NumberFormat.getIntegerInstance().format(t.size)}$$""" } else { "N/A" } ) }.toMap
    println(s"""${k} & ${methodsInOrder.map{ 
        mn => 
            var s = methodMap.get(mn).getOrElse("N/A") 
            if (mn == "mince" && s != "N/A") {
                s = s + " ($" + NumberFormat.getIntegerInstance().format(flipsMap(k).size) + "$) "
            }
            s
    }.mkString(" & ")}\\\\""")
    /*
    println(s"${k}\n==========")
    v.foreach{ case (fn, fi) => 
        println(s"""\t${fi.method} ${fi.size}""")
    }
    */
}

println("""\bottomrule""")
println("""\end{tabular}""")
println("""\end{center}""")


println()
import scala.util.parsing.json.JSONObject
import spray.json._
import DefaultJsonProtocol._ // !!! IMPORTANT, else `convertTo` and `toJson` won't work correctly

val minceResMap = minceMap.map{ case(k, v) => (k, v.size) }.toMap
val minceNoRCResMap = minceNoRCMap.map{ case(k, v) => (k, v.size) }.toMap
val fastqzResMap = fastqzMap.map{ case(k, v) => (k, v.size) }.toMap
val scalceResMap = scalceMap.map{ case(k, v) => (k, v.size) }.toMap
val flipsResMap = flipsMap.map{ case (k, v) => (k, v.size) }.toMap

//val jsonMap = Map("Mince (no RC)" -> minceNoRCResMap, "Mince" -> minceResMap).toJson //, "SCALCE" -> scalceResMap, "fastqz" -> fastqzResMap).toJson
val jsonMap = Map("flips" -> flipsResMap, "Mince (no RC)" -> minceNoRCResMap, "Mince" -> minceResMap,"SCALCE" -> scalceResMap, "fastqz" -> fastqzResMap).toJson
println(jsonMap)
