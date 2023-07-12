module.exports = {
  siteUrl: "https://cellxgene.cziscience.com/",
  generateRobotsTxt: true,
  exclude: ["/collections-sitemap.xml", "/cell-guide", "/wmg-v2"], // add additional directories or pages here if they should be excluded from the sitemap
  robotsTxtOptions: {
    policies: [
      // {userAgent: "*", disallow: "/folderToExclude/fileToExclude"} <-- use this option if you wish to exclude certain pages or folders from robots.txt
      {
        userAgent: "*",
        allow: "/",
        disallow: "https://api.cellxgene.cziscience.com/*",
      },
    ],
    additionalSitemaps: [
      "https://cellxgene.cziscience.com/collections-sitemap.xml",
    ],
  },
  siteMapSize: 50000,
};
