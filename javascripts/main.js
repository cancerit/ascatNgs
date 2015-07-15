var current_project = location.pathname.split("/")[1];
//var current_project = "cgpPindel"

var github_api_url = 'https://api.github.com/repos/cancerit/'.concat(current_project, '/readme?ref=master');

var corsproxy = "http://crossorigin.me/";
var github_url = "https://github.com/cancerit/".concat(current_project, "/raw/master/README.md");
var proxy_url = corsproxy.concat(github_url);

var error_html = "<small>Using fail-over page generation, please notify <a href='mailto:cgp-it@sanger.ac.uk' target='_top'>cgp-it@sanger.ac.uk</a></small><br>";

$.ajax({
	url:github_api_url,
	headers: {
	    "Accept": "application/vnd.github.3.html",
	    "User-Agent": "CancerIT"
	}
    })
    .success(function (response) {
	    $('#content').replaceWith(response); 
	})
    .error(function (r) {
	    $.get(proxy_url)
		.success(function (response) {
			var converter = new showdown.Converter(),
			    html = converter.makeHtml(response);
			$("#content").html(html);
		    });
	});

window.onload = function () {
    $("a.view_github_link" ).attr("href", "https://github.com/cancerit/".concat(current_project));
    $("a.zip_download_link").attr("href", "https://github.com/cancerit/".concat(current_project, "/zipball/master"));
    $("a.tar_download_link").attr("href", "https://github.com/cancerit/".concat(current_project, "/tarball/master"));
    $("h1.header").html(current_project);
    $("projectname").html(current_project);
    document.title = current_project;
};