<?php
$do = isset($_POST['do']) ? trim($_POST['do']) : '';
$id = isset($_POST['id']) ? intval($_POST['id']) : 0;
$job = isset($_POST['job']) ? trim($_POST['job']) : '';
if($do == 'search')
{
sleep(1);

$folder = "out/93u7pejy/";
if (strlen($job)>8){ exit;}
$folder = "out/".$job."/";
$nexttask ="";
	if(file_exists($folder."out.js"))
		{
			exit;
		}
	echo "<h4>You can bookmark this page and check your results later</h4> <h5>Job ID: $job</h5>";
	if(file_exists($folder."temp.geneout"))
		{
			echo "gene prediction done<br>";
			$nexttask = "conserved domain blast";
		}
	if(file_exists($folder."out.consDomain"))
		{
			echo "conserved domain blast done<br>";
			$nexttask ="finding homologs blast";
		}	

$k = 0;
$text= "";
while($k < $id){ $text .= "."; $k++;}
echo 'still waiting on '.$nexttask.$text." ".$id." seconds";
exit;
}
?>


<html>
<head>
<script type="text/javascript" src="js/jquery-1.7.2.min.js"></script>
<style type="text/css">
body {font-family:Arial, sans-serif}
</style>
</head>

<body>
<div id="result">Searching</div>
<script>
var id = 1;
var $_GET = getQueryParams(document.location.search);
$(document).ready(function()
{
	searchNext();
});
var job = $_GET['job'];

$('#again').on('click', function(e)
{
	e.preventDefault();
	id = 1;
	$('#result').html('Searching...');
	searchNext();
});

function getQueryParams(qs) {
    qs = qs.split("+").join(" ");
    var params = {},
        tokens,
        re = /[?&]?([^=]+)=([^&]*)/g;

    while (tokens = re.exec(qs)) {
        params[decodeURIComponent(tokens[1])]
            = decodeURIComponent(tokens[2]);
    }

    return params;
}

function searchNext()
{
$.ajax
({
type: 'POST',
url: '<?php echo $_SERVER['SCRIPT_NAME']; ?>',
data: 'do=search&id='+id+'&job='+job,
success: function(response)
{
if(!response)
{
$('<div id="end">Results found...redirecting now.   <META HTTP-EQUIV=Refresh CONTENT="1; URL=out/'+job+'/"> <a href="out/'+job+'">If you are not redirected, click here to manually fetch the page.</a>').insertAfter('#result');
}
else
{
id++;
$('#result').html(response);
return searchNext();
}
}
});
}
</script>
</body>
</html>
