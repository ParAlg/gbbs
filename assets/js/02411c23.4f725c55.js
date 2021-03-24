(window.webpackJsonp=window.webpackJsonp||[]).push([[3],{120:function(a,e,t){"use strict";t.d(e,"a",(function(){return i})),t.d(e,"b",(function(){return O}));var n=t(0),s=t.n(n);function m(a,e,t){return e in a?Object.defineProperty(a,e,{value:t,enumerable:!0,configurable:!0,writable:!0}):a[e]=t,a}function p(a,e){var t=Object.keys(a);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(a);e&&(n=n.filter((function(e){return Object.getOwnPropertyDescriptor(a,e).enumerable}))),t.push.apply(t,n)}return t}function r(a){for(var e=1;e<arguments.length;e++){var t=null!=arguments[e]?arguments[e]:{};e%2?p(Object(t),!0).forEach((function(e){m(a,e,t[e])})):Object.getOwnPropertyDescriptors?Object.defineProperties(a,Object.getOwnPropertyDescriptors(t)):p(Object(t)).forEach((function(e){Object.defineProperty(a,e,Object.getOwnPropertyDescriptor(t,e))}))}return a}function c(a,e){if(null==a)return{};var t,n,s=function(a,e){if(null==a)return{};var t,n,s={},m=Object.keys(a);for(n=0;n<m.length;n++)t=m[n],e.indexOf(t)>=0||(s[t]=a[t]);return s}(a,e);if(Object.getOwnPropertySymbols){var m=Object.getOwnPropertySymbols(a);for(n=0;n<m.length;n++)t=m[n],e.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(a,t)&&(s[t]=a[t])}return s}var b=s.a.createContext({}),l=function(a){var e=s.a.useContext(b),t=e;return a&&(t="function"==typeof a?a(e):r(r({},e),a)),t},i=function(a){var e=l(a.components);return s.a.createElement(b.Provider,{value:e},a.children)},N={inlineCode:"code",wrapper:function(a){var e=a.children;return s.a.createElement(s.a.Fragment,{},e)}},o=s.a.forwardRef((function(a,e){var t=a.components,n=a.mdxType,m=a.originalType,p=a.parentName,b=c(a,["components","mdxType","originalType","parentName"]),i=l(t),o=n,O=i["".concat(p,".").concat(o)]||i[o]||N[o]||m;return t?s.a.createElement(O,r(r({ref:e},b),{},{components:t})):s.a.createElement(O,r({ref:e},b))}));function O(a,e){var t=arguments,n=e&&e.mdxType;if("string"==typeof a||n){var m=t.length,p=new Array(m);p[0]=o;var r={};for(var c in e)hasOwnProperty.call(e,c)&&(r[c]=e[c]);r.originalType=a,r.mdxType="string"==typeof a?a:n,p[1]=r;for(var b=2;b<m;b++)p[b]=t[b];return s.a.createElement.apply(null,p)}return s.a.createElement.apply(null,t)}o.displayName="MDXCreateElement"},65:function(a,e,t){"use strict";t.r(e),t.d(e,"frontMatter",(function(){return p})),t.d(e,"metadata",(function(){return r})),t.d(e,"toc",(function(){return c})),t.d(e,"default",(function(){return l}));var n=t(3),s=t(7),m=(t(0),t(120)),p={id:"spanner",title:"Graph Spanner"},r={unversionedId:"benchmarks/sssp/spanner",id:"benchmarks/sssp/spanner",isDocsHomePage:!1,title:"Graph Spanner",description:"Problem Specification",source:"@site/docs/benchmarks/sssp/spanner.md",slug:"/benchmarks/sssp/spanner",permalink:"/gbbs/docs/benchmarks/sssp/spanner",version:"current",sidebar:"docs",previous:{title:"Single-Source Betweenness Centrality",permalink:"/gbbs/docs/benchmarks/sssp/ss_betweenness_centrality"},next:{title:"Biconnectivity",permalink:"/gbbs/docs/benchmarks/connectivity/biconnectivity"}},c=[{value:"Problem Specification",id:"problem-specification",children:[]},{value:"Algorithm Implementations",id:"algorithm-implementations",children:[]},{value:"Cost Bounds",id:"cost-bounds",children:[]},{value:"Compiling and Running",id:"compiling-and-running",children:[]},{value:"References",id:"references",children:[]}],b={toc:c};function l(a){var e=a.components,t=Object(s.a)(a,["components"]);return Object(m.b)("wrapper",Object(n.a)({},b,t,{components:e,mdxType:"MDXLayout"}),Object(m.b)("h2",{id:"problem-specification"},"Problem Specification"),Object(m.b)("h4",{id:"input"},"Input"),Object(m.b)("p",null,Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"G"),Object(m.b)("mo",{parentName:"mrow"},"="),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"V"),Object(m.b)("mo",{parentName:"mrow",separator:"true"},","),Object(m.b)("mi",{parentName:"mrow"},"E"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"G=(V, E)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"G"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(m.b)("span",{parentName:"span",className:"mrel"},"="),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"),Object(m.b)("span",{parentName:"span",className:"mpunct"},","),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.05764em"}},"E"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"))))),", an undirected, unweighted graph, and an integer\nstretch factor, ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"k")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"k")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.69444em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03148em"}},"k"))))),"."),Object(m.b)("h4",{id:"output"},"Output"),Object(m.b)("p",null,Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"H"),Object(m.b)("mo",{parentName:"mrow"},"\u2286"),Object(m.b)("mi",{parentName:"mrow"},"E")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"H \\subseteq E")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.8193em",verticalAlign:"-0.13597em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.08125em"}},"H"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(m.b)("span",{parentName:"span",className:"mrel"},"\u2286"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.05764em"}},"E"))))),", a set of edges such that for every ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"u"),Object(m.b)("mo",{parentName:"mrow",separator:"true"},","),Object(m.b)("mi",{parentName:"mrow"},"v"),Object(m.b)("mo",{parentName:"mrow"},"\u2208"),Object(m.b)("mi",{parentName:"mrow"},"V")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"u,v \\in V")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.7335400000000001em",verticalAlign:"-0.19444em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"u"),Object(m.b)("span",{parentName:"span",className:"mpunct"},","),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(m.b)("span",{parentName:"span",className:"mrel"},"\u2208"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.22222em"}},"V"))))),"\nconnected in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"G")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"G")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"0.68333em",verticalAlign:"0em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"G"))))),", ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("msub",{parentName:"mrow"},Object(m.b)("mrow",{parentName:"msub"},Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"d"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"i"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"s"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"t")),Object(m.b)("mi",{parentName:"msub"},"H")),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"u"),Object(m.b)("mo",{parentName:"mrow",separator:"true"},","),Object(m.b)("mi",{parentName:"mrow"},"v"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")"),Object(m.b)("mo",{parentName:"mrow"},"\u2264"),Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"k"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")"),Object(m.b)("mo",{parentName:"mrow"},"\u22c5"),Object(m.b)("msub",{parentName:"mrow"},Object(m.b)("mrow",{parentName:"msub"},Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"d"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"i"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"s"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"sans-serif"},"t")),Object(m.b)("mi",{parentName:"msub"},"G")),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"u"),Object(m.b)("mo",{parentName:"mrow",separator:"true"},","),Object(m.b)("mi",{parentName:"mrow"},"v"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"\\mathsf{dist}_{H}(u, v) \\leq O(k) \\cdot \\mathsf{dist}_{G}(u,v)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord"},Object(m.b)("span",{parentName:"span",className:"mord"},Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"d"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"i"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"s"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"t")),Object(m.b)("span",{parentName:"span",className:"msupsub"},Object(m.b)("span",{parentName:"span",className:"vlist-t vlist-t2"},Object(m.b)("span",{parentName:"span",className:"vlist-r"},Object(m.b)("span",{parentName:"span",className:"vlist",style:{height:"0.32833099999999993em"}},Object(m.b)("span",{parentName:"span",style:{top:"-2.5500000000000003em",marginRight:"0.05em"}},Object(m.b)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),Object(m.b)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},Object(m.b)("span",{parentName:"span",className:"mord mtight"},Object(m.b)("span",{parentName:"span",className:"mord mathdefault mtight",style:{marginRight:"0.08125em"}},"H"))))),Object(m.b)("span",{parentName:"span",className:"vlist-s"},"\u200b")),Object(m.b)("span",{parentName:"span",className:"vlist-r"},Object(m.b)("span",{parentName:"span",className:"vlist",style:{height:"0.15em"}},Object(m.b)("span",{parentName:"span"})))))),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"u"),Object(m.b)("span",{parentName:"span",className:"mpunct"},","),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),Object(m.b)("span",{parentName:"span",className:"mrel"},"\u2264"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03148em"}},"k"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),Object(m.b)("span",{parentName:"span",className:"mbin"},"\u22c5"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}})),Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord"},Object(m.b)("span",{parentName:"span",className:"mord"},Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"d"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"i"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"s"),Object(m.b)("span",{parentName:"span",className:"mord mathsf"},"t")),Object(m.b)("span",{parentName:"span",className:"msupsub"},Object(m.b)("span",{parentName:"span",className:"vlist-t vlist-t2"},Object(m.b)("span",{parentName:"span",className:"vlist-r"},Object(m.b)("span",{parentName:"span",className:"vlist",style:{height:"0.32833099999999993em"}},Object(m.b)("span",{parentName:"span",style:{top:"-2.5500000000000003em",marginRight:"0.05em"}},Object(m.b)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),Object(m.b)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},Object(m.b)("span",{parentName:"span",className:"mord mtight"},Object(m.b)("span",{parentName:"span",className:"mord mathdefault mtight"},"G"))))),Object(m.b)("span",{parentName:"span",className:"vlist-s"},"\u200b")),Object(m.b)("span",{parentName:"span",className:"vlist-r"},Object(m.b)("span",{parentName:"span",className:"vlist",style:{height:"0.15em"}},Object(m.b)("span",{parentName:"span"})))))),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"u"),Object(m.b)("span",{parentName:"span",className:"mpunct"},","),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03588em"}},"v"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"))))),"."),Object(m.b)("h2",{id:"algorithm-implementations"},"Algorithm Implementations"),Object(m.b)("p",null,"The algorithm is based on the Miller, Peng, Xu, and Vladu (MPXV) paper from SPAA'15 ","[1]",".\nOur implementation is described in more detail in ","[2]",".\nThe code for our implemenation is available\n",Object(m.b)("a",{parentName:"p",href:"https://github.com/ldhulipala/gbbs/tree/master/benchmarks/Spanner/MPXV15"},"here"),"."),Object(m.b)("p",null,"The construction results in an ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"k"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(k)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03148em"}},"k"),Object(m.b)("span",{parentName:"span",className:"mclose"},")"))))),"-spanner with expected size\n",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("msup",{parentName:"mrow"},Object(m.b)("mi",{parentName:"msup"},"n"),Object(m.b)("mrow",{parentName:"msup"},Object(m.b)("mn",{parentName:"mrow"},"1"),Object(m.b)("mo",{parentName:"mrow"},"+"),Object(m.b)("mn",{parentName:"mrow"},"1"),Object(m.b)("mi",{parentName:"mrow",mathvariant:"normal"},"/"),Object(m.b)("mi",{parentName:"mrow"},"k"))),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(n^{1+1/k})")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1.138em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord"},Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"n"),Object(m.b)("span",{parentName:"span",className:"msupsub"},Object(m.b)("span",{parentName:"span",className:"vlist-t"},Object(m.b)("span",{parentName:"span",className:"vlist-r"},Object(m.b)("span",{parentName:"span",className:"vlist",style:{height:"0.8879999999999999em"}},Object(m.b)("span",{parentName:"span",style:{top:"-3.063em",marginRight:"0.05em"}},Object(m.b)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),Object(m.b)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},Object(m.b)("span",{parentName:"span",className:"mord mtight"},Object(m.b)("span",{parentName:"span",className:"mord mtight"},"1"),Object(m.b)("span",{parentName:"span",className:"mbin mtight"},"+"),Object(m.b)("span",{parentName:"span",className:"mord mtight"},"1"),Object(m.b)("span",{parentName:"span",className:"mord mtight"},"/"),Object(m.b)("span",{parentName:"span",className:"mord mathdefault mtight",style:{marginRight:"0.03148em"}},"k"))))))))),Object(m.b)("span",{parentName:"span",className:"mclose"},")"))))),"."),Object(m.b)("h2",{id:"cost-bounds"},"Cost Bounds"),Object(m.b)("p",null,"The algorithm runs in ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"m"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(m)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"m"),Object(m.b)("span",{parentName:"span",className:"mclose"},")")))))," work and ",Object(m.b)("span",{parentName:"p",className:"math math-inline"},Object(m.b)("span",{parentName:"span",className:"katex"},Object(m.b)("span",{parentName:"span",className:"katex-mathml"},Object(m.b)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},Object(m.b)("semantics",{parentName:"math"},Object(m.b)("mrow",{parentName:"semantics"},Object(m.b)("mi",{parentName:"mrow"},"O"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},"("),Object(m.b)("mi",{parentName:"mrow"},"k"),Object(m.b)("mi",{parentName:"mrow"},"log"),Object(m.b)("mo",{parentName:"mrow"},"\u2061"),Object(m.b)("mi",{parentName:"mrow"},"n"),Object(m.b)("mo",{parentName:"mrow",stretchy:"false"},")")),Object(m.b)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"O(k\\log n)")))),Object(m.b)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},Object(m.b)("span",{parentName:"span",className:"base"},Object(m.b)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.02778em"}},"O"),Object(m.b)("span",{parentName:"span",className:"mopen"},"("),Object(m.b)("span",{parentName:"span",className:"mord mathdefault",style:{marginRight:"0.03148em"}},"k"),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mop"},"lo",Object(m.b)("span",{parentName:"span",style:{marginRight:"0.01389em"}},"g")),Object(m.b)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),Object(m.b)("span",{parentName:"span",className:"mord mathdefault"},"n"),Object(m.b)("span",{parentName:"span",className:"mclose"},")")))))," depth w.h.p.\nMore details about our implementation and the cost bounds can be found in Section 6.1 of ","[2]","."),Object(m.b)("h2",{id:"compiling-and-running"},"Compiling and Running"),Object(m.b)("p",null,"The benchmark can be compiled by running:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"bazel build -c opt //benchmarks/Spanner/MPXV15:Spanner\n")),Object(m.b)("p",null,"It can then be run on a test input graph in the ",Object(m.b)("em",{parentName:"p"},"uncompressed format")," as follows:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"numactl -i all ./bazel-bin/benchmarks/Spanner/MPXV15/Spanner_main -s -m -src 1 inputs/rMatGraph_J_5_100\n")),Object(m.b)("p",null,"It can then be run on a test input graph in the ",Object(m.b)("em",{parentName:"p"},"compressed format")," as follows:"),Object(m.b)("pre",null,Object(m.b)("code",{parentName:"pre"},"numactl -i all ./bazel-bin/benchmarks/Spanner/MPXV15/Spanner_main -s -c -m -src 1 inputs/rMatGraph_J_5_100.bytepda\n")),Object(m.b)("h2",{id:"references"},"References"),Object(m.b)("p",null,"[1]"," Gary L. Miller, Richard Peng, Adrian Vladu, Shen Chen Xu.",Object(m.b)("br",null),"\n",Object(m.b)("a",{parentName:"p",href:"https://arxiv.org/abs/1309.3545"},Object(m.b)("em",{parentName:"a"},"Improved Parallel Algorithms for Spanners and Hopsets")),". ",Object(m.b)("br",null),"\nProceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), 2015."),Object(m.b)("p",null,"[2]"," Laxman Dhulipala, Guy Blelloch, and Julian Shun",Object(m.b)("br",null),"\n",Object(m.b)("a",{parentName:"p",href:"https://ldhulipala.github.io/papers/gbbs_topc.pdf"},Object(m.b)("em",{parentName:"a"},"Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable")),Object(m.b)("br",null),"\nProceedings of the ACM Symposium on Parallelism in Algorithms and Architectures (SPAA), pp. 393-404, 2018. ",Object(m.b)("br",null)))}l.isMDXComponent=!0}}]);