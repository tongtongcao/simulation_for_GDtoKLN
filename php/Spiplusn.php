<?php
	echo "\n";
	for ($i = 0; $i < 86; $i++) {
		$url_left = 'http://gwdac.phys.gwu.edu/cgi-bin/go3pr2?sl=CM12&rt=2&ot=S&iv=A&il=0&ii=3&iu=180&fv=E&fn=';
		$url_right = '1000&jpeg=PLOT&u=995&l=1005';
		$E = 900+20*$i;
		$url = $url_left.substr_replace($url_right, strval($E), 0, 4);
		echo "\n";
		echo "Fetching S for E=";
		echo $E;
		echo " ...";
		echo "\n";
		echo $url;
		echo "\n";
		$content = file_get_contents($url, NULL, NULL, 0, 1000000);
		$outfilename = '0000S.txt';
		file_put_contents(substr_replace($outfilename, strval($E), 0, 4), substr($content, 0, 1000000));
	}
	echo "\n";
?>

