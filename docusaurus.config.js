const remarkMath = require('remark-math');
const rehypeKatex = require('rehype-katex');

module.exports = {
  title: 'GBBS',
  tagline: 'Graph Based Benchmark Suite',
  url: 'https://paralg.github.io',
  favicon: 'img/graph.svg',
  organizationName: 'ParAlg', // Usually your GitHub org/user name.
  projectName: 'gbbs', // Usually your repo name.
  url: 'https://paralg.github.io',
  baseUrl: '/gbbs/',
  themeConfig: {
   //  disableDarkMode: true,
    navbar: {
      title: 'GBBS',
/*      style: 'light', */
      logo: {
        alt: 'GBBS Logo',
        src: 'img/graph.svg',
      },
      links: [
        {
          to: 'docs/introduction',
          activeBasePath: 'docs',
          label: 'Docs',
          position: 'left',
        },
        {
          to: 'docs/library/overview',
          activeBasePath: 'docs',
          label: 'Library',
          position: 'left',
        },
        {
          to: 'docs/benchmarks/overview',
          activeBasePath: 'docs/',
          label: 'Benchmark Implementations',
          position: 'left',
        },
        {
          to: 'https://github.com/ParAlg/gbbs',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    footer: {
/*      style: 'light', */
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Tutorial',
              to: 'docs/tutorial',
            },
            {
              label: 'Documentation',
              to: 'docs/introduction',
            },
            {
              label: 'Library',
              to: 'docs/library/overview',
            },
            {
              label: 'Benchmark Specifications',
              to: 'docs/benchmarks/overview',
            }
          ],
        },
        {
          title: 'Community',
          items: [
            {
              label: 'Contributors',
              to: 'docs/contributors',
            },
            {
              label: 'Publications and Resources',
              to: 'docs/research',
            },
          ],
        },
        {
          title: 'Links',
          items: [
            {
              label: 'GitHub',
              to: 'https://github.com/ParAlg/gbbs',
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} GBBS Authors.`,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          remarkPlugins: [remarkMath],
          rehypePlugins: [rehypeKatex],
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
  stylesheets: [
    {
      href: 'https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.css',
      type: 'text/css',
    }
  ],
/*  scripts: [
    {
      src: 'https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.js',
      async: true,
    },
  ], */
};

//<script defer src="https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.js" integrity="sha384-y23I5Q6l+B6vatafAwxRu/0oK/79VlbSz7Q9aiSZUvyWYIYsd+qj+o24G5ZU2zJz" crossorigin="anonymous"></script>
//    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.11.1/dist/katex.min.css" integrity="sha384-zB1R0rpPzHqg7Kpt0Aljp8JPLqbXI3bhnPWROx27a9N0Ll6ZP/+DiW/UqRcLbRjq" crossorigin="anonymous">
