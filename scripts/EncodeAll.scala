#!/bin/bash
exec scala -nc -savecompiled "$0" "$@"
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


val retMince = getFileList(new File(args(0)), "mince", List("offs.lz", "seqs.lz"))

val minceMap = new HMap[String, FileInfo]
retMince.foreach{ case (k,v) => 
    var base = k.getName.split("""\.""")(0)
    if (minceMap.contains(base)) {
        minceMap(base).size += v.size
    } else {
        minceMap.put(base, v)
    }
}

val scalceMap = new HMap[String, FileInfo]
val fastqzMap = new HMap[String, FileInfo]
/*
val retScalce = getFileList(new File(args(1)), "scalce", "scalcer")

retScalce.foreach{ case (k,v) => 
    var base = k.getName.split("_")(0)
    if (scalceMap.contains(base)) {
        scalceMap(base).size += v.size
        scalceMap(base).isPE = true
    } else {
        scalceMap.put(base, v)
    }
}

val retFastqz = getFileList(new File(args(2)), "fastqz", "fxb.zpaq")

retFastqz.foreach{ case (k,v) =>
    var base = k.getName.split("""\.""")(0)
    if (fastqzMap.contains(base)) {
        fastqzMap(base).size += v.size
    } else {
        fastqzMap.put(base, v)
    }
}
*/
val mergedMap = (minceMap.toSeq ++ scalceMap.toSeq ++ fastqzMap.toSeq).groupBy{ case (k,v) => k}

val numKeys = mergedMap.keySet.size
val methodsInOrder = Array("fastqz", "scalce", "mince")

println("""\begin{center}""")
println("""\begin{tabular}{lrrr}""")

println("""\toprule""")
println(s"""Dataset & ${methodsInOrder.mkString(" & ")}\\\\""")
println("""\midrule""")

mergedMap.foreach{ case (k,v) =>
    val methodMap = v.map{ case (n,t) => 
        (t.method, if (t.size > 0) { s"""$$${NumberFormat.getIntegerInstance().format(t.size)}$$""" } else { "N/A" } ) }.toMap
    println(s"""${k} & ${methodsInOrder.map{ mn => methodMap.get(mn).getOrElse("N/A") }.mkString(" & ")}\\\\""")
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
