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
            <title>Privacy - CZ CELLxGENE</title>
          </Head>
          <header>
            <Image
              data-testid="cellxgene-logo"
              src={rawCellxgeneLogo}
              alt="CELLxGENE logo"
              width="119"
              height="35"
            />
          </header>

          <br />

          <div>
            <h1>Privacy Policy</h1>
            <p>Last updated: July 25, 2023.</p>
            {/* Introduction */}
            <>
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
                Please read our full{" "}
                <a href={ROUTES.TOS}>Terms of Use (&quot;Terms&quot;)</a> (which
                incorporates this Privacy Policy) and this Privacy Policy for
                complete details, but here is the key information you should
                know:
              </p>
              <ul className="text-list">
                <li>
                  This is a <strong>free service</strong> we provide in order to
                  advance biomedical research.
                </li>
                <li>
                  The datasets made available through the Site are{" "}
                  <strong>not personally identifiable</strong>.
                </li>
                <li>
                  We use the privacy-friendly Plausible service to collect{" "}
                  <strong>basic analytics</strong> about our Site traffic (e.g.
                  the number of visitors and page views) so we know how it’s
                  being used. We also analyze the usage of the Site which helps
                  us improve the Site as well as get a sense of its impact.
                </li>
              </ul>
            </>
            {/* Data Controllers and Contracting Parties */}
            <>
              <h2>Data Controllers and Contracting Parties</h2>
              <p>
                By accessing and using the Site, you are contracting with the
                Chan Zuckerberg Initiative Foundation, a 501(c)(3) nonprofit
                private foundation (“<strong>Provider</strong>,” “
                <strong>we</strong>” or “<strong>us</strong>”), and agreeing
                that Provider is the “controller” of your personal data provided
                to, collected by, or processed in connection with the Site. This
                Privacy Policy along with the{" "}
                <a href={ROUTES.TOS}>Terms of Use</a> form a contract.{" "}
                <strong>
                  If you don’t agree with this Privacy Policy or the Terms of
                  Use, do not access or use the Site.
                </strong>{" "}
                This Privacy Policy applies to only this Site, and excludes any
                other services that state that they are offered under a
                different privacy policy. For example, this Privacy Policy does
                not apply to chanzuckerberg.com.
              </p>
            </>
            {/* Data We Collect */}
            <>
              <h2 id="data-we-collect">Data We Collect</h2>
              <p>
                Cellxgene is a tool that enables fast visualizations of curated
                single-cell transcriptomics datasets. It takes data submitted by
                researchers (<em>gene expression matrices</em> that indicate
                counts of how often genes are expressed in certain cell types
                and <em>metadata</em> detailing how those matrices were
                generated) and helps you visualize it for faster analysis,
                exploration, and – hopefully – insight.
              </p>
              <p>
                The datasets available on the Site are not personally
                identifiable.
              </p>
              <ul className="text-list">
                <li>
                  <h4>Data You Provide To Us.</h4> We collect certain
                  information from you which may include identifiers or
                  professional or employment-related information. This includes
                  data you provide us as part of a submission or registering for
                  an account (ex: name and email address). This information is
                  necessary in order to create your account and provide you with
                  access to the services we offer on the Site. We collect your
                  information if you decide to sign up for our newsletter, if
                  you contact us for support or with information about a dataset
                  via email, or if you respond to a survey or participate in a
                  feedback session. We also collect the data you submit for
                  display in CELLxGENE (ex: single-cell .h5ad matrix files).
                </li>
                <li>
                  <h4>Data From Your Browser or Device.</h4> Whenever you use
                  any online service, certain internet or other electronic
                  network activity information gets created and logged
                  automatically; the same is true when you access or use the
                  Site. Here’s what we collect:
                  <ul>
                    <li>
                      <h5>Log.</h5> When you access or use the Site (whether on
                      your computer or on a mobile device), we gather certain
                      information automatically and store it in log files. This
                      information includes IP addresses, the Internet Site
                      Provider, referring pages, date/time stamps, clickstream
                      data, and duration of time spent on the Site. We collect
                      this data pursuant to our legitimate interests in
                      understanding how Visitors use our Site, improving our
                      Site, and securing our Site.
                    </li>
                    <li>
                      <h5>Device.</h5> In addition to log data, we collect
                      information about the device you’re using to access the
                      Site; this includes the type of device, browser type, and
                      operating system that helps us understand when something
                      goes wrong. We collect this data pursuant to our
                      legitimate interests in improving our Site.
                    </li>
                    <li>
                      <h5>Cookies and Other Similar Technologies.</h5> For
                      certain users that need to log into CELLxGENE Discover to
                      publish data, we use essential cookies (small text files
                      sent by your computer each time you access the Site that
                      are unique to your account or your browser) to enable that
                      use. We also use local storage to save your preferences
                      (e.g., to remember that you blocked a pop-up). For web
                      analytics, we do not use Google Analytics. Instead, we use
                      the privacy-friendly Plausible as our website analytic
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
                  </ul>
                </li>
                <li>
                  <h4>Other Sources.</h4> We may display publicly available
                  information about you when you mention us on social media
                  networks, such as Twitter. This is reflected in the “CELL X
                  GENE IN THE NEWS” section on the Site.
                </li>
                <li>
                  <h4>Sensitive Personal Information.</h4> We do not collect
                  sensitive personal information about you, and thus we do not
                  use or disclose your sensitive personal information.
                </li>
              </ul>
            </>
            {/* How We Use Your Data */}
            <>
              <h2>How We Use Your Data</h2>
              <p>
                CZIF does not sell your personal information nor do we share
                your personal data for behavioral advertising purposes. We do
                use your data for the following business purposes:
              </p>
              <ul className="text-list">
                <li>
                  <h4>Site.</h4> We use the information we collect to provide,
                  maintain, secure, and improve the Site, including
                  understanding the content that Visitors find valuable.
                </li>
                <li>
                  <h4>Communications.</h4> We may also use your information to
                  respond to an email from you, send a survey, engage about a
                  dataset you wish to share, or send our newsletter if you
                  opt-in to receiving it. Survey responses and feedback sessions
                  help us understand and improve the Site.
                </li>
                <li>
                  <h4>Aggregate Insights and Analytics.</h4> We use the
                  privacy-friendly Plausible service to produce and share
                  aggregated insights that do not identify you. For example we
                  may receive from Plausible statistics about the location of
                  our Visitors, and how many Visitors engage with the Site on a
                  monthly basis. These aggregated insights are not personally
                  identifiable. We also use your information to analyze how our
                  Site is being used, such as understanding how many downloads
                  there have been for a specific dataset.
                </li>
                <li>
                  <h4>Other Sources.</h4> We use information from Other Sources
                  to understand what you find valuable, what we can improve, and
                  what impact our Site is having on the science community.
                </li>
                <li>
                  <h4>Security and Investigations.</h4> We use your data
                  (including your communications) if we think it’s necessary for
                  security purposes or to investigate violations of our{" "}
                  <a href={ROUTES.TOS}>Terms of Use</a> or this Privacy Policy.
                </li>
              </ul>
            </>
            {/* Retention and Deletion */}
            <>
              <h2>Retention and Deletion</h2>
              <p>
                We retain your personal data as long as we believe we need it
                for the purposes we have collected it or to meet legal
                obligations, resolve disputes, maintain security, prevent fraud
                and abuse, or enforce our agreements with you. This includes
                data you provided to us and data generated or inferred from your
                use of the Site.
              </p>
            </>
            {/* How We Disclose Information */}
            <>
              <h2>How We Disclose Information</h2>
              <p>
                Except in the instances listed below, we will not disclose your
                personal information to others unless you consent to it, nor
                will we ever sell your personal information to advertisers or
                other third parties. However, we disclose your information in
                the following ways:
              </p>
              <ul className="text-list">
                <li>
                  <h4>Third Party Site Providers.</h4> CZIF works with service
                  providers that help us operate, secure, and improve the Site.
                  These services are, for example, performing statistical
                  analysis, database management services, database hosting,
                  survey providers, and security. To the extent they will have
                  access to your information, their use is limited by this
                  Privacy Policy.
                </li>
                <li>
                  <h4>Legal and Safety Reasons.</h4> We may disclose information
                  if we believe in good faith that it’s necessary (a) in
                  connection with any legal investigation; (b) to comply with
                  relevant laws or to respond to subpoenas or warrants served on
                  us; (c) to protect or defend our rights or property; and/or
                  (d) to investigate or assist in preventing any violation of
                  the law.
                </li>
                <li>
                  <h4>CZIF Entities and Affiliates.</h4> The Chan Zuckerberg
                  Initiative, LLC (“CZI LLC”) is our primary technology partner,
                  focusing on the Site’s infrastructure, security, and
                  compliance. In this role, CZI LLC is a data controller for all
                  data referenced in this Privacy Policy. As with Service
                  Providers mentioned above, CZI LLC&lsquo;s use of data is
                  limited by this Privacy Policy. “Affiliates” refers to
                  entities controlled by or under common control with CZIF (such
                  as CZI LLC) and does not include Meta Platforms, Inc. for
                  purposes of this policy.
                </li>
                <li>
                  <h4>Reorganization, Sale or Merger.</h4> We may disclose your
                  information in connection with a merger, reorganization, or
                  sale of all or a portion of our organization or assets related
                  to CZIF. In the event of a merger, reorganization or sale of
                  assets, the buyer or other successor entity will continue to
                  be bound by the terms of this Privacy Policy.
                </li>
              </ul>
            </>
            {/* Choices and Rights */}
            <>
              <h2>Choices and Rights</h2>
              <div className="inline">
                <h4>Rights.</h4> You have the following rights with respect to
                the personal data we have about you:
              </div>
              <ul>
                <li>
                  <h5>Delete data.</h5> You can ask us to erase or delete all or
                  some of your personal data subject to our legal obligations
                  and lawful exceptions.
                </li>
                <li>
                  <h5>Change or correct personal data.</h5> You can also ask us
                  to change, update, or fix inaccurate data in certain cases,
                  subject to our legal obligations and lawful exceptions.
                </li>
                <li>
                  <h5>Object to, limit, or restrict use of personal data.</h5>{" "}
                  You can ask us to stop using all or some of your personal data
                  (e.g., if we have no legal right to keep using it) or to limit
                  our use of it (e.g., if your personal data is inaccurate or
                  unlawfully held).
                </li>
                <li>
                  <h5>Right to access and/or take your personal data.</h5> You
                  can ask us for a copy of your personal data in
                  machine-readable form.
                </li>
                <li>
                  <h5>Right to notice.</h5> You have a right to receive notice
                  of our personal information collection, use, retention, and
                  disclosure practices at or before collection of personal
                  information.
                </li>
                <li>
                  <h5>The right not to be discriminated against.</h5> CZIF will
                  not discriminate against you in any manner for exercising any
                  of the above rights with respect to your personal data.
                </li>
              </ul>
              <p>
                If you would like to exercise your right to any of the above,
                email us at{" "}
                <a href="mailto:privacy@chanzuckerberg.com">
                  privacy@chanzuckerberg.com
                </a>
                . In the email, please provide us with your name, the country
                (and state if within the United States) in which you live, which
                of the above rights you would like to exercise, and sufficient
                information that allows us to reasonably verify that you are the
                person about whom we collected personal information. If you
                would like an authorized agent to make a request for you, have
                that agent email{" "}
                <a href="mailto:privacy@chanzuckerberg.com">
                  privacy@chanzuckerberg.com
                </a>{" "}
                with the above information along with additional information
                sufficient for us to verify that the authorized agent is acting
                on your behalf. Please also let us know if you have questions or
                concerns related to exercising any rights you have under
                applicable law to control your personal data.
              </p>
              <p>
                If you would like to appeal a CZI decision with respect to a
                request to exercise any of these rights, please email us at{" "}
                <a href="mailto:privacy@chanzuckerberg.com">
                  privacy@chanzuckerberg.com
                </a>{" "}
                and explain the basis for your appeal.
              </p>
              <p>
                If you wish to raise a concern about our use of your information
                (and without prejudice to any other rights you may have), you
                have the right to do so with your local supervisory authority.
              </p>
            </>
            {/* Security */}
            <>
              <h2>Security</h2>
              <p>
                Security of personal data is important to us. We implement
                security safeguards designed to protect your personal data,
                including reasonable administrative, technical, and physical
                safeguards to protect that personal data from unauthorized
                access, use, alteration and destruction. Despite these efforts,
                we cannot guarantee that your data may not be accessed,
                disclosed, altered, or destroyed by a breach of any of our
                physical, technical, or administrative safeguards. Please notify
                us immediately at{" "}
                <a href="mailto:security@chanzuckerberg.com">
                  security@chanzuckerberg.com
                </a>{" "}
                if you become aware of any security issues relating to the Site.
              </p>
            </>
            {/* Data Transfers */}
            <>
              <h2>Data Transfers</h2>
              <p>
                CZIF is based in the United States; when you engage with the
                Site, you are sending personal data into the United States which
                may have different data protection rules than those of your
                country. We process data both inside and outside of the United
                States.
              </p>
            </>
            {/* Our Legal Bases */}
            <>
              <h2>Our Legal Bases</h2>
              <p>
                We will collect, use and disclose your personal data only where
                we have a legal right to do so. This section explains our legal
                bases for processing personal data, including under GDPR.
              </p>
              <ul className="text-list">
                <li>
                  <h4>Consent.</h4> We rely on consent to engage in certain data
                  collection activities, like data you choose to submit to us
                  for display on CELLxGENE Discover or if you decide to sign up
                  for our newsletter.
                </li>
                <li>
                  <h4>Legitimate Interests.</h4> We rely on legitimate interests
                  to process the device and log data we collect when you use the
                  Site. We process this data based on our legitimate interest in
                  understanding how the Site is being used so we can improve it
                  and have a sense of its impact, and your legitimate interest
                  in accessing the Site.
                </li>
                <li>
                  <h4>Contract.</h4> We rely on contract where processing is
                  necessary for the performance of a contract with you (e.g., to
                  make the data you submitted to us publicly available via
                  CELLxGENE Discover).
                  <p>
                    Where we rely on consent, you have the right to revoke your
                    consent and where we rely on legitimate interests, you have
                    the right to object by emailing us at{" "}
                    <a href="mailto:privacy@chanzuckerberg.com">
                      privacy@chanzuckerberg.com
                    </a>
                    . If you have any questions about the lawful bases on which
                    we collect and use your personal data, please contact us at{" "}
                    <a href="mailto:GDPR-REP@chanzuckerberg.com">
                      GDPR-REP@chanzuckerberg.com
                    </a>
                    .
                  </p>
                </li>
              </ul>
            </>
            {/* Additional Information For California Residents */}
            <>
              <h2>Additional Information For California Residents</h2>
              <p>
                The California Consumer Privacy Act (“CCPA”) requires certain
                businesses to give California residents a number of rights
                regarding their personal information. We are offering these
                rights to you, including the right to have your personal
                information deleted (subject to certain exceptions), the right
                to change or correct your personal information, the right to
                limit the use or disclosure of your sensitive personal
                information (if applicable), the right to access your personal
                information, the right to opt-out of the “selling” or “sharing”
                of personal information (if applicable), and the right not to be
                discriminated against for exercising these rights.
              </p>
              <p>
                These rights, and how to exercise them, are described in more
                detail in the Section titled “Choices and Rights” of this Policy
                Policy. In addition to these rights, we give you a right to
                request the following information about your personal
                information that we have collected in the past 12 months:
              </p>
              <>
                <div className="inline">
                  <h4>The Right to Know.</h4> This right allows you to request
                  the following information about the personal information that
                  we’ve collected about you in the past 12 months:
                </div>
                <ul>
                  <li>
                    <h5>Information about Data Collection</h5>
                    <ul>
                      <li>
                        The categories of personal information that have been
                        collected about you.
                      </li>
                      <li>
                        The categories of sources from which we have collected
                        personal information.
                      </li>
                      <li>
                        The business purpose for which we have collected
                        personal information.
                      </li>
                    </ul>
                  </li>
                  <li>
                    <h5>Information about Data Disclosure</h5>
                    <ul>
                      <li>
                        The categories of personal information, if any, that
                        have been sold, shared, or disclosed for a business
                        purpose to third parties.
                      </li>
                      <li>
                        The categories of third parties to whom personal
                        information was sold, shared, or disclosed for a
                        business purpose.
                      </li>
                      <li>
                        Identification of the specific business purpose for
                        disclosing the consumer’s personal information.
                      </li>
                    </ul>
                  </li>
                </ul>
              </>
              <>
                <p>
                  We have described in fuller detail in this Privacy Policy the
                  personal information that we collect, how we use, and disclose
                  it, but provide the following additional disclosure:
                </p>
                <h4>Information about Data Collection</h4>
                <ul>
                  <li>
                    <h5>Information we collect.</h5> We have collected the
                    following categories of personal information from consumers
                    within the past 12 months: (1) identifiers; (2) professional
                    or employment-related information; (3) internet or other
                    electronic network activity; (4) geolocation data; and (5)
                    information provided through survey responses or feedback
                    sessions. We have not sold this information nor have we
                    shared this information for behavioral advertising purposes.
                  </li>
                  <li>
                    <h5>Sources of information.</h5> We obtain these categories
                    of personal information from the sources described in the
                    Section titled “Data We Collect.”
                  </li>
                  <li>
                    <h5>Purposes of collection.</h5> We collect personal
                    information for one or more of the following business
                    purposes as described in the Section titled “How We Use Your
                    Data” above.
                  </li>
                </ul>
                <h4>Information about Data Disclosure</h4>
                <ul>
                  <li>
                    <h5>Information we disclose.</h5> We have disclosed the
                    following categories of personal information within the past
                    12 months: (1) identifiers; (2) professional or
                    employment-related information; (3) internet or other
                    electronic network activity; (4) geolocation data; and (5)
                    information provided through survey responses or feedback
                    sessions. We have not sold this information nor have we
                    shared this information for behavioral advertising purposes.
                  </li>
                  <li>
                    <h5>Third parties to whom we disclose.</h5> The categories
                    of third parties to whom we have disclosed this personal
                    information are described in the Section titled “How We
                    Disclose Information” in this Privacy Policy.
                  </li>
                  <li>
                    <h5>Purposes of disclosure.</h5> We disclose the personal
                    information we collect about you for one or more of the
                    following business purposes described in the Section titled
                    “How We Use Your Data” above.
                  </li>
                  <li>
                    <h5>Sensitive data.</h5> We do not collect sensitive
                    personal information, and thus we do not disclose or use
                    your sensitive information.
                  </li>
                </ul>
              </>
            </>
            {/* Additional Information for Residents of Virginia, Colorado,
                Connecticut, and Utah */}
            <>
              <h2>
                Additional Information for Residents of Virginia, Colorado,
                Connecticut, and Utah
              </h2>
              <p>
                Virginia, Colorado, Connecticut, and Utah also have adopted
                privacy laws that give consumers certain rights, including the
                right to confirm whether controllers are processing the
                consumer’s personal data, the right to access that data, the
                right to obtain a copy of that data, the right to correct
                inaccuracies in that data, and the right to delete that data. As
                discussed above in the Section titled “Choices and Rights,” we
                provide these rights to all consumers, regardless of where they
                reside.
              </p>
              <p>
                Additionally, these four states have adopted rights to opt-out
                of: (1) targeted advertising; (2) the sale of personal data; and
                (3) profiling in furtherance of decisions that produce legal or
                similarly significant effects concerning the consumer. We do not
                sell your data, use it for targeted advertising or to profile
                you in furtherance of decisions that produce legal or similarly
                significant effects.
              </p>
            </>
            {/* Other Important Information */}
            <>
              <h2>Other Important Information</h2>
              <ul className="text-list">
                <li>
                  <h4>Direct Marketing and Do Not Track Signals.</h4> We don’t
                  currently disclose personal data with third parties for their
                  direct marketing purposes, nor do we support any Do Not Track
                  (&quot;DNT&quot;) signals, since there’s currently no standard
                  for how online services respond to those signals. As standards
                  develop, we may establish policies for responding to DNT
                  signals that we would describe in this Privacy Policy.
                </li>
                <li>
                  <h4>Children.</h4> We do not have actual knowledge that we
                  have sold or shared the personal information of users under 16
                  years of age.
                </li>
                <li>
                  <h4>Changes.</h4> We may modify this Privacy Policy from time
                  to time, and you can see when the last update was by looking
                  at the “Last Updated” date at the top of this Policy. If we
                  make material changes to it, we’ll provide you notice through
                  this Privacy Policy. Your continued use of the Site after we
                  publish a notice about changes to this Privacy Policy means
                  that you acknowledge and agree to the updated Privacy Policy
                  following the date it takes effect.
                </li>
                <li>
                  <h4>Contact Information.</h4> If you have questions or
                  complaints regarding this Privacy Policy, please contact us at{" "}
                  <a href="mailto:privacy@chanzuckerberg.com">
                    privacy@chanzuckerberg.com
                  </a>
                  .
                  <div>
                    To comply with article 27 of the GDPR and the UK-GDPR, we
                    have appointed a representative who can accept
                    communications on behalf of CZIF and CZI LLC in relation to
                    personal data processing activities falling within the scope
                    of the GDPR or the UK-GDPR. If you wish to contact them,
                    their details are as follows:
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
                </li>
              </ul>
            </>
          </div>
        </PrivacyStyle>
      </CommonStyle>
    </Layout>
  );
};

export default Privacy;
