import React from 'react';
import classnames from 'classnames';
import Layout from '@theme/Layout';
import Link from '@docusaurus/Link';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';
import useBaseUrl from '@docusaurus/useBaseUrl';
import styles from './styles.module.css';

const features = [
  {
    title: <>Functional Graph Interface</>,
    imageUrl: 'img/graph-figure.svg',
    description: (
      <>
        The GBBS library provides a graph-processing interface that
        extends the Ligra, Ligra+, and Julienne graph processing
        interfaces with other useful functional primitives, including
        map, reduce, and filter, defined over both vertices and
        graphs. The API also supports compression completely under the
        hood.
      </>
    ),
  },
  {
    title: <>Fundamental Graph Problems</>,
    imageUrl: 'img/undraw_docusaurus_tree.svg',
    description: (
      <>
        The GBBS benchmark, implemented using the GBBS library,
        provides provably-efficient implementations for over 22
        fundamental graph problems including shortest-path problems,
        connectivity problems, and many data-mining problems. Good
        theoretical bounds ensure consistent performance across
        different inputs.
      </>
    ),
  },
  {
    title: <>State-of-the-Art Performance</>,
    imageUrl: 'img/undraw_docusaurus_react.svg',
    description: (
      <>
        GBBS implementations achieve state-of-the-art performance. For
        example, our implementation of connectivity solves the
        WebDataCommons Hyperlink2012 graph, the largest publicly-available graph with
        over 3 billion vertices and over 200 billion edges in 8.5
        seconds on a 72-core machine.
      </>
    ),
  },
];

function Feature({imageUrl, title, description}) {
  const imgUrl = useBaseUrl(imageUrl);
  return (
    <div className={classnames('col col--4', styles.feature)}>
      {imgUrl && (
        <div className="text--center">
          <img className={styles.featureImage} src={imgUrl} alt={title} />
        </div>
      )}
      <h3>{title}</h3>
      <p>{description}</p>
    </div>
  );
}

function Home() {
  const context = useDocusaurusContext();
  const {siteConfig = {}} = context;
  return (
    <Layout
      title={`${siteConfig.title}`}
      description="Description will go into a meta tag in <head />">
      <header className={classnames('hero', styles.heroBanner)}>
        <div className="container">
          <h1 className="hero__title">{siteConfig.title}</h1>
          <p className="hero__subtitle">{siteConfig.tagline}</p>
          <div className={styles.buttons}>
            <Link
              className={classnames(
                'button button--outline button--secondary button--lg',
                styles.getStarted,
              )}
              to={useBaseUrl('docs/introduction')}>
              Get Started
            </Link>
          </div>
        </div>
      </header>
      <main>
        {features && features.length && (
          <section className={styles.features}>
            <div className="container">
              <div className="row">
                {features.map((props, idx) => (
                  <Feature key={idx} {...props} />
                ))}
              </div>
            </div>
          </section>
        )}
      </main>
    </Layout>
  );
}

export default Home;
