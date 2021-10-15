import Head from "next/head";
import Image from "next/image";
import { ROUTES } from "src/common/constants/routes";
import rawCellxgeneLogo from "src/components/common/staticPages/cellxgene.png";
import {
  CommonStyle,
  Layout,
  PrivacyStyle,
} from "src/components/common/staticPages/style";

const Privacy = (): JSX.Element => {
  return (
    <Layout>
      <CommonStyle>
        <PrivacyStyle>
          <Head>
            <title>cellxgene | Privacy</title>
          </Head>
          <header>
            <Image
              data-test-id="cellxgene-logo"
              src={rawCellxgeneLogo}
              alt="cellxgene logo"
              width="119"
              height="35"
            />
          </header>

          <br />

          <main>
            <h1>Privacy Policy</h1>

            <p>Last updated: Oct 19, 2021.</p>

            <br />

            <b>Summary of key changes in this Oct 19, 2021 update:</b>

            <br />

            <span>
              We no longer use Google Analytics for web analytics. Instead, we
              have transitioned to use the privacy-friendly Plausible service to
              collect basic analytics about our Site traffic to understand how
              it’s used. As a result, we collect less personal data about you.
            </span>

            <br />
            <br />

            <h2>Introduction</h2>

            <p>
              We are a site for biologists, computational researchers, and
              others (collectively, “<strong>Visitors</strong>” or “
              <strong>you</strong>
              ”) to explore curated single-cell transcriptomics datasets at
              cellxgene.cziscience.com (the “<strong>Site</strong>”). As a
              visitor to the Site, the collection, use and sharing of your
              personal data is subject to this Privacy Policy.
            </p>

            <p>
              Please read our full <a href={ROUTES.TOS}>Terms</a> (which
              incorporates this Privacy Policy) and this Privacy Policy for
              complete details, but here is the key information you should know:
            </p>

            <ul>
              <li>
                This is a <strong>free service</strong> we provide in order to
                advance biomedical research.
              </li>
              <li>
                The data made available through the Site is{" "}
                <strong>not personally identifiable</strong>.
              </li>
              <li>
                We use the privacy-friendly Plausible service to collect{" "}
                <strong>basic analytics</strong> about our Site traffic (e.g.
                the number of visitors and page views) so we know how it’s being
                used. This helps us improve the Site as well as get a sense of
                its impact.
              </li>
            </ul>

            <ol>
              <li>
                <h3>Data Controllers and Contracting Parties</h3>
                <p>
                  By accessing and using the Site, you are contracting with the
                  Chan Zuckerberg Initiative Foundation, a 501(c)(3) nonprofit
                  private foundation (“<strong>Provider</strong>,” “
                  <strong>we</strong>” or “<strong>us</strong>”), and agreeing
                  that Provider is the “controller” of your personal data
                  provided to, collected by, or processed in connection with the
                  Site. This Privacy Policy along with the{" "}
                  <a href={ROUTES.TOS}>Terms of Use</a> form a contract.{" "}
                  <strong>
                    If you don’t agree with this Privacy Policy or the Terms of
                    Use, do not access or use the Site.
                  </strong>{" "}
                  This Privacy Policy applies to only this Site, and excludes
                  any other services that state that they are offered under a
                  different privacy policy. For example, this Privacy Policy
                  does not apply to chanzuckerberg.com.
                </p>
              </li>

              <li>
                <h3>Data We Collect</h3>
                <p>
                  Cellxgene is a tool that enables fast visualizations of
                  curated single-cell transcriptomics datasets. It takes data
                  submitted by researchers (<em>gene expression matrices</em>{" "}
                  that indicate counts of how often genes are expressed in
                  certain cell types and <em>metadata</em> detailing how those
                  matrices were generated) and helps you visualize it for faster
                  analysis, exploration, and – hopefully – insight.
                </p>
                <p>This data is not personally identifiable.</p>
                <ul className="section2">
                  <p>
                    <h4>Data You Provide To Us.</h4> We collect certain
                    information from you when you provide it to us directly.
                    Specifically, this includes data you submit for display in
                    cellxgene (ex: single-cell .h5ad matrix files), data you
                    provide you us as part of submission or registering for an
                    account (ex: name and email address), as well as if you
                    contact us for support or with information about a dataset
                    via email. This information is necessary in order to create
                    your account and provide you with access to the services we
                    offer on the Site.
                  </p>
                  <p>
                    <h4>Data From Your Browser or Device.</h4> Whenever you use
                    any online service, certain information gets created and
                    logged automatically; the same is true when you access or
                    use the Site. Here’s what we collect:
                    <ul>
                      <p>
                        <li>
                          <h5>Log.</h5> When you access or use the Site (whether
                          on your computer or on a mobile device), we gather
                          certain information automatically and store it in log
                          files. This information includes IP addresses, the
                          Internet Site Provider, referring pages, date/time
                          stamps, clickstream data, and duration of time spent
                          on the Site. We collect this data pursuant to our
                          legitimate interests in understanding how Visitors use
                          our Site and improving our Site.
                        </li>
                        <li>
                          <h5>Device.</h5> In addition to log data, we collect
                          information about the device you’re using to access
                          the Site; this includes the type of device, browser
                          type, and operating system that helps us understand
                          when something goes wrong. We collect this data
                          pursuant to our legitimate interests in improving our
                          Site.
                        </li>
                        <li>
                          <h5>Cookies and Other Similar Technologies.</h5> For
                          certain users that need to log into the cellxgene
                          portal to publish data, we use essential cookies
                          (small text files sent by your computer each time you
                          access the Site that are unique to your account or
                          your browser) to enable that use. For web analytics,
                          we do not use Google Analytics. Instead, we use the
                          privacy-friendly Plausible as our website analytic
                          tool. Learn more about Plausible’s data and privacy
                          practices{" "}
                          <a
                            href="https://plausible.io/data-policy"
                            rel="noopener"
                            target="_blank"
                          >
                            here
                          </a>
                          .
                        </li>
                      </p>
                    </ul>
                  </p>
                </ul>
              </li>
              <li>
                <h3>How We Use Your Data</h3>

                <p>
                  CZI does not sell, rent, or lease your personal data to
                  others. We use your data for the following business purposes
                  only:
                </p>

                <ol className="section3">
                  <p>
                    <h4>Site.</h4> We use the information we collect to provide
                    the Site, and maintain and improve the Site, including
                    understanding the content that Visitors find valuable.
                  </p>
                  <p>
                    <h4>Communications.</h4> We may also use your information to
                    respond to an email from you, and to engage with you about a
                    dataset you wish to share.
                  </p>
                  <p>
                    <h4>Aggregate Insights.</h4> We use the privacy-friendly
                    Plausible service to produce and share aggregated insights
                    that do not identify you. For example we may receive from
                    Plausible statistics about the location of our Visitors, and
                    how many Visitors engage with the Site on a monthly basis.
                    These aggregated insights are not personally identifiable.
                  </p>
                  <p>
                    <h4>Security and Investigations.</h4> We use your data
                    (including your communications) if we think it’s necessary
                    for security purposes or to investigate violations of our{" "}
                    <a href={ROUTES.TOS}>Terms of Use</a> or this Privacy
                    Policy. We may use human and automated systems and
                    inferences we make to determine whether you or others can be
                    trusted to engage with the Sites.
                  </p>
                </ol>
              </li>

              <li>
                <h3>How We Share Information</h3>

                <p>
                  Except in the instances listed below, we will not disclose
                  your personal information to others unless you consent to it,
                  nor will we ever sell your personal information to advertisers
                  or other third parties. However, we share your information in
                  the following ways:
                </p>

                <ol className="section4">
                  <p>
                    <h4>Third Party Site Providers.</h4> CZIF works with service
                    providers that help us operate, secure, and improve the
                    Site. These services are, for example, performing
                    statistical analysis, database management services, database
                    hosting, and security. To the extent they will have access
                    to your information, their use is limited by this Privacy
                    Policy.
                  </p>
                  <p>
                    <h4>Legal and Safety Reasons.</h4> We may disclose
                    information if we believe in good faith that it’s necessary
                    (a) in connection with any legal investigation; (b) to
                    comply with relevant laws or to respond to subpoenas or
                    warrants served on us; (c) to protect or defend our rights
                    or property; and/or (d) to investigate or assist in
                    preventing any violation of the law.
                  </p>
                  <p>
                    <h4>CZIF Entities and Affiliates.</h4> The Chan Zuckerberg
                    Initiative, LLC (“CZI LLC”) is our primary technology
                    partner, focusing on the Site’s infrastructure, security,
                    and compliance. In this role, CZI LLC is a data controller
                    for all data referenced in this Privacy Policy. As with
                    Service Providers mentioned above, CZI LLC&lsquo;s use of
                    data is limited by this Privacy Policy. “Affiliates” refers
                    to entities controlled by or under common control with CZIF
                    (such as CZI LLC) and does not include Facebook for purposes
                    of this policy.
                  </p>
                  <p>
                    <h4>Reorganization, Sale or Merger.</h4> We may share your
                    information in connection with a merger, reorganization, or
                    sale of all or a portion of our organization or assets
                    related to CZIF. In the event of a merger, reorganization or
                    sale of assets, the buyer or other successor entity will
                    continue to be bound by the terms of this Privacy Policy.
                  </p>
                </ol>
              </li>
              <li>
                <h3>Choices and Rights</h3>
                <ol className="section5">
                  <p>
                    <h4>Rights.</h4> You have the following rights with respect
                    to the personal data we have about you:
                  </p>
                  <p>
                    <ul>
                      <li>
                        <h5>Delete data.</h5> You can ask us to erase or delete
                        all or some of your personal data.
                      </li>
                      <li>
                        <h5>Change or correct personal data.</h5> You can also
                        ask us to change, update or fix your data in certain
                        cases, particularly if it’s inaccurate.
                      </li>
                      <li>
                        <h5>
                          Object to, limit, or restrict use of personal data.
                        </h5>{" "}
                        You can ask us to stop using all or some of your
                        personal data (e.g., if we have no legal right to keep
                        using it) or to limit our use of it (e.g., if your
                        personal data is inaccurate or unlawfully held).
                      </li>
                      <li>
                        <h5>Right to access and/or take your personal data.</h5>{" "}
                        You can ask us for a copy of your personal data in
                        machine-readable form.
                      </li>
                      <li>
                        <h5>The right not to be discriminated against.</h5> CZIF
                        will not discriminate against you in any manner for
                        exercising any of the above rights with respect to your
                        personal data.
                      </li>
                      <li>
                        Contact us at{" "}
                        <a href="mailto:privacy@chanzuckerberg.com">
                          privacy@chanzuckerberg.com
                        </a>{" "}
                        if you have questions or would like to exercise any
                        rights you have under applicable law to control your
                        personal data. If you wish to raise a concern about our
                        use of your information (and without prejudice to any
                        other rights you may have), you have the right to do so
                        with your local supervisory authority.
                      </li>
                    </ul>
                  </p>
                </ol>
              </li>

              <li>
                <h3>Retention and Deletion.</h3>
                <p>
                  We retain your personal data as needed to provide the Site.
                  This includes data you provided to us and data generated or
                  inferred from your use of the Site.
                </p>
              </li>

              <li>
                <h3>Data Transfers.</h3>
                <p>
                  CZIF is based in the United States; when you engage with the
                  Site, you are sending personal data into the United States
                  which may have different data protection rules than those of
                  your country. We process data both inside and outside of the
                  United States.
                </p>
              </li>

              <li>
                <h3>Our Legal Bases.</h3>
                <p>
                  We will collect, use and share your personal data only where
                  we have a legal right to do so. This section explains our
                  legal bases for processing personal data, including under
                  GDPR.
                </p>
                <ol className="section8">
                  <p>
                    <h4>Consent.</h4> We rely on consent to engage in certain
                    data collection activities, like data you choose to submit
                    to us for display on cellxgene.
                  </p>
                  <p>
                    <h4>Legitimate Interests.</h4> We rely on legitimate
                    interests to process the device and log data we collect when
                    you use the Site. We process this data based on our
                    legitimate interest in understanding how the Site is being
                    used so we can improve it and have a sense of its impact,
                    and your legitimate interest in accessing the Site.
                  </p>
                  <p>
                    <h4>Contract.</h4> We rely on contract where processing is
                    necessary for the performance of a contract with you (e.g.
                    to make the data you submitted to us publicly available via
                    cellxgene).
                  </p>
                </ol>
                <p>
                  Where we rely on consent, you have the right to revoke your
                  consent and where we rely on legitimate interests, you have
                  the right to object by emailing us at{" "}
                  <a href="mailto:privacy@chanzuckerberg.com">
                    privacy@chanzuckerberg.com
                  </a>
                  . If you have any questions about the lawful bases on which we
                  collect and use your personal data, please contact our Data
                  Protection Officer via email at{" "}
                  <a href="mailto:GDPR-REP@chanzuckerberg.com">
                    GDPR-REP@chanzuckerberg.com
                  </a>
                  .
                </p>
              </li>

              <li>
                <h3>Other Important Information</h3>
                <ol className="section9">
                  <p>
                    <h4>Security.</h4> Security of personal data is important to
                    us. We implement security safeguards designed to protect
                    your personal data, including reasonable administrative,
                    technical and physical safeguards to protect that personal
                    data from unauthorized access, use, alteration and
                    destruction. Despite these efforts, we cannot guarantee that
                    your data may not be accessed, disclosed, altered, or
                    destroyed by a breach of any of our physical, technical, or
                    administrative safeguards. Please notify us immediately at{" "}
                    <a href="mailto:security@chanzuckerberg.com">
                      security@chanzuckerberg.com
                    </a>{" "}
                    if you become aware of any security issues relating to the
                    Site.
                  </p>
                  <p>
                    <h4>Direct Marketing and Do Not Track Signals.</h4> We don’t
                    currently share personal data with third parties for their
                    direct marketing purposes, nor do we support any Do Not
                    Track signals, since there’s currently no standard for how
                    online services respond to those signals. As standards
                    develop, we may establish policies for responding to DNT
                    signals that we would describe in this Privacy Policy.
                  </p>
                  <p>
                    <h4>Changes.</h4> We may modify this Privacy Policy from
                    time to time, and you can see when the last update was by
                    looking at the “Last Updated” date at the top of this page.
                    If we make material changes to it, we’ll provide you notice
                    through this Privacy Policy. Your continued use of the Site
                    after we publish a notice about changes to this Privacy
                    Policy means that you acknowledge and agree to the updated
                    Privacy Policy following the date it takes effect.
                  </p>
                  <p>
                    <h4>Children.</h4> The Site is not designed or intended for
                    children under 16. If we become aware that we have the
                    information of such children collected through the Site, we
                    will promptly delete it.
                  </p>
                  <p>
                    <h4>Contact Information.</h4> If you have questions or
                    complaints regarding this Privacy Policy, please contact us
                    at{" "}
                    <a href="mailto:privacy@chanzuckerberg.com">
                      privacy@chanzuckerberg.com
                    </a>
                    .
                    <div>
                      To comply with article 27 of the GDPR and the UK-GDPR, we
                      have appointed a representative who can accept
                      communications on behalf of CZIF and CZI LLC in relation
                      to personal data processing activities falling within the
                      scope of the GDPR or the UK-GDPR. If you wish to contact
                      them, their details are as follows:
                      <br />
                      <br />
                      <address>
                        European GDPR Representative: <br />
                        Bird & Bird GDPR Representative Services SRL <br />
                        Avenue Louise 235 <br />
                        1050 Bruxelles <br />
                        Belgium
                        <br />
                        <a href="mailto:EUrepresentative.ChanZuckerberg@twobirds.com">
                          EUrepresentative.ChanZuckerberg@twobirds.com
                        </a>
                      </address>
                      <address>
                        UK Data Protection Representative:
                        <br />
                        Bird & Bird GDPR Representative Services UK
                        <br />
                        12 New Fetter Lane <br />
                        London EC4A 1JP <br />
                        United Kingdom
                        <br />
                        <a href="mailto:UKrepresentative.ChanZuckerberg@twobirds.com">
                          UKrepresentative.ChanZuckerberg@twobirds.com
                        </a>
                      </address>
                    </div>
                  </p>
                </ol>
              </li>
            </ol>
          </main>
        </PrivacyStyle>
      </CommonStyle>
    </Layout>
  );
};

export default Privacy;
