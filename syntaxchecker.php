<?
function base36_encode($base10){
    return base_convert($base10,10,36);
}

$now = time();
$random = rand(0,8);
$random2 = rand(2,9);
//print "now ".$random2.$now.$random."<br>";
$token =  base36_encode($random2.$now.$random);



function run_in_background($Command, $Priority = 0)
   {
       if($Priority)
           $PID = shell_exec("nohup nice -n $Priority $Command > /dev/null 2>&1 &");
       else
           $PID = shell_exec("$Command > /dev/null 2>&1 &");
       return($PID);
   }
   function is_process_running($PID)
   {
       exec("ps $PID", $ProcessState);
       return(count($ProcessState) >= 2);
   }

function flush_buffers(){ 
    ob_end_flush(); 
    ob_flush(); 
    flush(); 
    ob_start(); 
} 


$filename = "Inputs/".$token.".fa";
$fh = fopen($filename, 'w');
$data = $_POST["inputSeq"];
fwrite($fh, $data);
fclose($fh);




$src = "HTML_include/";
$dest = "out/$token/";

mkdir($dest);

$include1 = "index.html";
$include2 = "raphael-min.js";
$include3 = "json2.js";
copy($src.$include1, $dest.$include1);
copy($src.$include2, $dest.$include2);
copy($src.$include3, $dest.$include3);

print "<META HTTP-EQUIV=Refresh CONTENT=\"1; URL=status.php?job=$token\">";


$command="python SyntaxChecker.py $filename out/$token > out/$token/errlog 2>&1 &";
$ps = run_in_background($command);


exit;
?>
