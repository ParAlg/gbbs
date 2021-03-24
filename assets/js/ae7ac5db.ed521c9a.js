(window.webpackJsonp=window.webpackJsonp||[]).push([[31],{101:function(e,n,t){"use strict";t.r(n),t.d(n,"frontMatter",(function(){return o})),t.d(n,"metadata",(function(){return i})),t.d(n,"toc",(function(){return p})),t.d(n,"default",(function(){return s}));var r=t(3),a=(t(0),t(120));const o={id:"python_bindings",title:"Python Bindings"},i={unversionedId:"python_bindings",id:"python_bindings",isDocsHomePage:!1,title:"Python Bindings",description:"Next, we illustrate how to use GBBS from Python using the python",source:"@site/docs/python_bindings.md",slug:"/python_bindings",permalink:"/gbbs/docs/python_bindings",version:"current",sidebar:"docs",previous:{title:"Graph Formats",permalink:"/gbbs/docs/formats"},next:{title:"BFS",permalink:"/gbbs/docs/tutorial/bfs_tutorial"}},p=[{value:"Loading SNAP Graphs",id:"loading-snap-graphs",children:[]},{value:"Loading Graphs in the Adjacency Array Format",id:"loading-graphs-in-the-adjacency-array-format",children:[]},{value:"Running Graph Algorithms",id:"running-graph-algorithms",children:[]}],c={toc:p};function s({components:e,...n}){return Object(a.b)("wrapper",Object(r.a)({},c,n,{components:e,mdxType:"MDXLayout"}),Object(a.b)("p",null,"Next, we illustrate how to use GBBS from Python using the python\nbindings we provide (through ",Object(a.b)("a",{parentName:"p",href:"https://pybind11.readthedocs.io/en/stable/"},"pybind11"),")"),Object(a.b)("p",null,"Extending the Python bindings after implementing a new benchmark\nrequires only a few lines of code to add an extra method to the graph\nobject exported by the library.\nNote that the contents of this file have been tested using Python\n3.7.6. You may run into issues using an earlier version."),Object(a.b)("p",null,"We first build the bindings using Bazel and add the compiled libraries\nto the Python path:"),Object(a.b)("pre",null,Object(a.b)("code",{parentName:"pre",className:"language-sh"},"$ bazel build //pybindings:gbbs_lib.so\n$ export PYTHONPATH=$(pwd)/bazel-bin/pybindings/:$PYTHONPATH\n")),Object(a.b)("h2",{id:"loading-snap-graphs"},"Loading SNAP Graphs"),Object(a.b)("p",null,"Launch the Python REPL, import the library, and import a downloaded\ngraph from the SNAP dataset."),Object(a.b)("pre",null,Object(a.b)("code",{parentName:"pre",className:"language-python"},'>>> import gbbs\n>>> G = gbbs.loadSnap("com-youtube.ungraph.txt", undirected=True)\n')),Object(a.b)("p",null,"This command creates an uncompressed graph in the GBBS format at the\nsame location as the input (compression can optionally be enabled\nusing a separate flag)."),Object(a.b)("h2",{id:"loading-graphs-in-the-adjacency-array-format"},"Loading Graphs in the Adjacency Array Format"),Object(a.b)("p",null,"We also provide commands to load graphs that have already been created\nin some valid format."),Object(a.b)("p",null,"To load a graph in text format:"),Object(a.b)("pre",null,Object(a.b)("code",{parentName:"pre",className:"language-python"},'>>> import gbbs\n>>> G = gbbs.loadGraph("rMatGraph_J_5_100", undirected=True, compressed=False, binary=True)\n')),Object(a.b)("p",null,"To load a graph in compressed format:"),Object(a.b)("pre",null,Object(a.b)("code",{parentName:"pre",className:"language-python"},'>>> import gbbs\n>>> G = gbbs.loadGraph("rMatGraph_J_5_100.bytepda", undirected=True, compressed=True, binary=False)\n')),Object(a.b)("p",null,"To load a graph in the binary format:"),Object(a.b)("pre",null,Object(a.b)("code",{parentName:"pre",className:"language-python"},'>>> import gbbs\n>>> G = gbbs.loadGraph("rMatGraph_J_5_100.binary", undirected=True, compressed=False, binary=True)\n')),Object(a.b)("h2",{id:"running-graph-algorithms"},"Running Graph Algorithms"),Object(a.b)("p",null,"We can then apply the methods defined on graphs as follows:"),Object(a.b)("pre",null,Object(a.b)("code",{parentName:"pre"},">>> sims = G.PageRank()\n")),Object(a.b)("p",null,"Other primitives can be applied similarly. For example:"),Object(a.b)("pre",null,Object(a.b)("code",{parentName:"pre"},">>> components = G.Connectivity()\n>>> print(components[10] == components[82])\nTrue\n>>> cores = G.KCore() # Computes coreness values\n>>> print(cores[10], cores[82])\n(41, 50)\n")))}s.isMDXComponent=!0},120:function(e,n,t){"use strict";t.d(n,"a",(function(){return b})),t.d(n,"b",(function(){return h}));var r=t(0),a=t.n(r);function o(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function i(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);n&&(r=r.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,r)}return t}function p(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?i(Object(t),!0).forEach((function(n){o(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):i(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function c(e,n){if(null==e)return{};var t,r,a=function(e,n){if(null==e)return{};var t,r,a={},o=Object.keys(e);for(r=0;r<o.length;r++)t=o[r],n.indexOf(t)>=0||(a[t]=e[t]);return a}(e,n);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(r=0;r<o.length;r++)t=o[r],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(a[t]=e[t])}return a}var s=a.a.createContext({}),l=function(e){var n=a.a.useContext(s),t=n;return e&&(t="function"==typeof e?e(n):p(p({},n),e)),t},b=function(e){var n=l(e.components);return a.a.createElement(s.Provider,{value:n},e.children)},d={inlineCode:"code",wrapper:function(e){var n=e.children;return a.a.createElement(a.a.Fragment,{},n)}},u=a.a.forwardRef((function(e,n){var t=e.components,r=e.mdxType,o=e.originalType,i=e.parentName,s=c(e,["components","mdxType","originalType","parentName"]),b=l(t),u=r,h=b["".concat(i,".").concat(u)]||b[u]||d[u]||o;return t?a.a.createElement(h,p(p({ref:n},s),{},{components:t})):a.a.createElement(h,p({ref:n},s))}));function h(e,n){var t=arguments,r=n&&n.mdxType;if("string"==typeof e||r){var o=t.length,i=new Array(o);i[0]=u;var p={};for(var c in n)hasOwnProperty.call(n,c)&&(p[c]=n[c]);p.originalType=e,p.mdxType="string"==typeof e?e:r,i[1]=p;for(var s=2;s<o;s++)i[s]=t[s];return a.a.createElement.apply(null,i)}return a.a.createElement.apply(null,t)}u.displayName="MDXCreateElement"}}]);