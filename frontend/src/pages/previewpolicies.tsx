import React, { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import rawCellxgeneLogo from "src/components/common/staticPages/cellxgene.png";
import rawCZILogo from "src/components/common/staticPages/CZI_Logotype_RGB.png";
import {
  CellxgeneLogo,
  CommonStyle,
  CZILogo,
  Layout,
  PrivacyStyle,
} from "src/components/common/staticPages/style";
import SEO from "src/components/seo";

const PreviewPolicies: FC = () => {
  return (
    <Layout>
      <CommonStyle>
        <PrivacyStyle>
          <SEO title="Preview Policies" />
          <header>
            <CZILogo
              data-test-id="czi-logo"
              src={String(rawCZILogo)}
              alt="CZI logo"
            />
            <CellxgeneLogo
              data-test-id="cellxgene-logo"
              src={String(rawCellxgeneLogo)}
              alt="cellxgene logo"
            />
          </header>

          <h1>
            Summary of April 2021 Updates to Terms of Service and Privacy Policy
          </h1>

          <p>
            <span>
              We periodically update our terms and policies to ensure that they
              are transparent and easy to understand. Here is a summary of the
              key updates from April 2021:
              <ul>
                <li>
                  Cellxgene is now provided by the Chan Zuckerberg Initiative
                  Foundation. The Chan Zuckerberg Initiative has moved its
                  Science team to its 501(c)(3) private foundation from CZI,
                  LLC. We’re committed to preserving the same experience for
                  you. For questions about these changes, please contact us at{" "}
                  <a href="mailto:privacy@chanzuckerberg.com">
                    privacy@chanzuckerberg.com
                  </a>
                  .
                </li>
              </ul>
            </span>
          </p>
        </PrivacyStyle>
      </CommonStyle>

      <Layout>
        <CommonStyle>
          <PrivacyStyle>
            {/* <SEO title="Privacy" /> */}
            <header>
              <CZILogo
                data-test-id="czi-logo"
                src={String(rawCZILogo)}
                alt="CZI logo"
              />
              <CellxgeneLogo
                data-test-id="cellxgene-logo"
                src={String(rawCellxgeneLogo)}
                alt="cellxgene logo"
              />
            </header>

            <br />
            <h3>
              Jump to <a href="#tos">Terms of Use</a>.
            </h3>
            <br />

            <main>
              <h1>Privacy Policy</h1>

              <p>Last updated: April 1, 2021.</p>

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
                complete details, but here is the key information you should
                know:
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
                  We collect <strong>basic analytics</strong> about your usage
                  of the Site (e.g. frequency, duration, IP address) so we know
                  how it’s being used. This helps us improve the Site as well as
                  get a sense of its impact.
                </li>
              </ul>

              <ol>
                <li>
                  <h3>Data Controllers and Contracting Parties</h3>
                  <p>
                    By accessing and using the Site, you are contracting with
                    the Chan Zuckerberg Initiative Foundation, a 501(c)(3)
                    nonprofit private foundation (“<strong>Provider</strong>,” “
                    <strong>we</strong>” or “<strong>us</strong>”), and agreeing
                    that Provider is the “controller” of your personal data
                    provided to, collected by, or processed in connection with
                    the Site. This Privacy Policy along with the{" "}
                    <a href={ROUTES.TOS}>Terms of Use</a> form a contract.{" "}
                    <strong>
                      If you don’t agree with this Privacy Policy or the Terms
                      of Use, do not access or use the Site.
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
                    matrices were generated) and helps you visualize it for
                    faster analysis, exploration, and – hopefully – insight.
                  </p>
                  <p>This data is not personally identifiable.</p>
                  <ol className="section2">
                    <li>
                      <h4>Data You Provide To Us.</h4> We collect certain
                      information from you when you provide it to us directly.
                      Specifically, this includes data you submit for display in
                      cellxgene (ex: single-cell .loom matrix files), data you
                      provide you us as part of submission or registering for an
                      account (ex: name and email address), as well as if you
                      contact us for support or with information about a dataset
                      via email. This information is necessary in order to
                      create your account and provide you with access to the
                      services we offer on the Site.
                    </li>
                    <li>
                      <h4>Data From Your Browser or Device.</h4> Whenever you
                      use any online service, certain information gets created
                      and logged automatically; the same is true when you access
                      or use the Site. Here’s what we collect:
                      <ul>
                        <li>
                          <h5>Log.</h5> When you access or use the Site (whether
                          on your computer or on a mobile device), we gather
                          certain information automatically and store it in log
                          files. This information includes IP addresses, the
                          Internet Site Provider, referring pages, date/time
                          stamps, clickstream data, and duration of time spent
                          on the Site. We collect this data pursuant to our
                          legitimate interests in understanding how Visitors use
                          our Site.
                        </li>
                        <li>
                          <h5>Device.</h5> In addition to log data, we collect
                          information about the device you’re using to access
                          the Site; this includes the type of device, browser
                          type, operating system, settings, unique device
                          identifiers, and crash data that helps us understand
                          when something goes wrong. We collect this data
                          pursuant to our legitimate interests in improving our
                          Site.
                        </li>
                        <li>
                          <h5>Cookies and Other Similar Technologies.</h5> We
                          also use cookies (small text files sent by your
                          computer each time you access the Site that are unique
                          to your account or your browser) and similar
                          technologies. For example, we use analytics services
                          (e.g.{" "}
                          <a href="https://www.google.com/url?q=https://www.google.com/policies/privacy/partners/&amp;sa=D&amp;ust=1585348324239000">
                            Google Analytics
                          </a>
                          ) that place cookies that collect information that
                          allows us to understand how often you use the Sites,
                          where you are accessing the Sites from and events that
                          happen on the Sites. If you’d like to opt-out of
                          Google Analytics, you can go{" "}
                          <a href="https://www.google.com/url?q=https://tools.google.com/dlpage/gaoptout&amp;sa=D&amp;ust=1585348324240000">
                            here
                          </a>
                          .
                        </li>
                      </ul>
                    </li>
                  </ol>
                </li>
                <li>
                  <h3>How We Use Your Data</h3>

                  <p>
                    CZI does not sell, rent, or lease your personal data to
                    others. We use your data for the following business purposes
                    only:
                  </p>

                  <ol className="section3">
                    <li>
                      <h4>Site.</h4> We use the information we collect to
                      provide the Site, and maintain and improve the Site,
                      including understanding the content that Visitors find
                      valuable.
                    </li>
                    <li>
                      <h4>Communications.</h4> We may also use your information
                      to respond to an email from you, and to engage with you
                      about a dataset you wish to share.
                    </li>
                    <li>
                      <h4>Aggregate Insights.</h4> We use your data to produce
                      and share aggregated insights that do not identify you.
                      For example we may use your data to generate statistics
                      about the location of our Visitors, and how many Visitors
                      engage with the Site on a monthly basis. These aggregated
                      insights are not personally identifiable.
                    </li>
                    <li>
                      <h4>Security and Investigations.</h4> We use your data
                      (including your communications) if we think it’s necessary
                      for security purposes or to investigate violations of our{" "}
                      <a href={ROUTES.TOS}>Terms of Use</a> or this Privacy
                      Policy. We may use human and automated systems and
                      inferences we make to determine whether you or others can
                      be trusted to engage with the Sites.
                    </li>
                  </ol>
                </li>

                <li>
                  <h3>How We Share Information</h3>

                  <p>
                    Except in the instances listed below, we will not disclose
                    your personal information to others unless you consent to
                    it, nor will we ever sell your personal information to
                    advertisers or other third parties. However, we share your
                    information in the following ways:
                  </p>

                  <ol className="section4">
                    <li>
                      <h4>Third Party Site Providers.</h4> CZIF works with
                      service providers that help us operate, secure, and
                      improve the Site. These services are, for example,
                      performing statistical analysis, database management
                      services, database hosting, and security. To the extent
                      they will have access to your information, their use is
                      limited by this Privacy Policy.
                    </li>
                    <li>
                      <h4>Legal and Safety Reasons.</h4> We may disclose
                      information if we believe in good faith that it’s
                      necessary (a) in connection with any legal investigation;
                      (b) to comply with relevant laws or to respond to
                      subpoenas or warrants served on us; (c) to protect or
                      defend our rights or property; and/or (d) to investigate
                      or assist in preventing any violation of the law.
                    </li>
                    <li>
                      <h4>CZIF Entities and Affiliates.</h4> The Chan Zuckerberg
                      Initiative, LLC (“CZI LLC”) is our primary technology
                      partner, focusing on the Site’s infrastructure, security,
                      and compliance. In this role, CZI LLC is a data controller
                      for all data referenced in this Privacy Policy. As with
                      Service Providers mentioned above, CZI LLC&lsquo;s use of
                      data is limited by this Privacy Policy. “Affiliates”
                      refers to entities controlled by or under common control
                      with CZIF (such as CZI LLC) and does not include Facebook
                      for purposes of this policy.
                    </li>
                    <li>
                      <h4>Reorganization, Sale or Merger.</h4> We may share your
                      information in connection with a merger, reorganization,
                      or sale of all or a portion of our organization or assets
                      related to CZIF. In the event of a merger, reorganization
                      or sale of assets, the buyer or other successor entity
                      will continue to be bound by the terms of this Privacy
                      Policy.
                    </li>
                  </ol>
                </li>
                <li>
                  <h3>Choices and Rights</h3>
                  <ol className="section5">
                    <li>
                      <h4>Rights.</h4> You have the following rights with
                      respect to the personal data we have about you:
                    </li>
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
                    <li>
                      <h4>Consent.</h4> We rely on consent to engage in certain
                      data collection activities, like through cookies or data
                      you choose to submit to us for display on cellxgene.
                    </li>
                    <li>
                      <h4>Legitimate Interests.</h4> We rely on legitimate
                      interests to process the device and log data we collect
                      when you use the Site. We process this data based on our
                      legitimate interest in understanding how the Site is being
                      used so we can improve it and have a sense of its impact,
                      and your legitimate interest in accessing the Site.
                    </li>
                    <li>
                      <h4>Contract.</h4> We rely on contract where processing is
                      necessary for the performance of a contract with you (e.g.
                      to make the data you submitted to us publicly available
                      via cellxgene).
                    </li>
                  </ol>
                  <p>
                    Where we rely on consent, you have the right to revoke your
                    consent and where we rely on legitimate interests, you have
                    the right to object by emailing us at{" "}
                    <a href="mailto:privacy@chanzuckerberg.com">
                      privacy@chanzuckerberg.com
                    </a>
                    . If you have any questions about the lawful bases on which
                    we collect and use your personal data, please contact our
                    Data Protection Officer via email at{" "}
                    <a href="mailto:GDPR-REP@chanzuckerberg.com">
                      GDPR-REP@chanzuckerberg.com
                    </a>
                    .
                  </p>
                </li>

                <li>
                  <h3>Other Important Information</h3>
                  <ol className="section9">
                    <li>
                      <h4>Security.</h4> Security of personal data is important
                      to us. We implement security safeguards designed to
                      protect your personal data, including reasonable
                      administrative, technical and physical safeguards to
                      protect that personal data from unauthorized access, use,
                      alteration and destruction. Despite these efforts, we
                      cannot guarantee that your data may not be accessed,
                      disclosed, altered, or destroyed by a breach of any of our
                      physical, technical, or administrative safeguards. Please
                      notify us immediately at{" "}
                      <a href="mailto:security@chanzuckerberg.com">
                        security@chanzuckerberg.com
                      </a>{" "}
                      if you become aware of any security issues relating to the
                      Site.
                    </li>
                    <li>
                      <h4>Direct Marketing and Do Not Track Signals.</h4> We
                      don’t currently share personal data with third parties for
                      their direct marketing purposes, nor do we support any Do
                      Not Track signals, since there’s currently no standard for
                      how online services respond to those signals. As standards
                      develop, we may establish policies for responding to DNT
                      signals that we would describe in this Privacy Policy.
                    </li>
                    <li>
                      <h4>Changes.</h4> We may modify this Privacy Policy from
                      time to time, and you can see when the last update was by
                      looking at the “Last Updated” date at the top of this
                      page. If we make material changes to it, we’ll provide you
                      notice through this Privacy Policy. Your continued use of
                      the Site after we publish a notice about changes to this
                      Privacy Policy means that you acknowledge and agree to the
                      updated Privacy Policy following the date it takes effect.
                    </li>
                    <li>
                      <h4>Children.</h4> The Site is not designed or intended
                      for children under 16. If we become aware that we have the
                      information of such children collected through the Site,
                      we will promptly delete it.
                    </li>
                    <li>
                      <h4>Contact Information.</h4> If you have questions or
                      complaints regarding this Privacy Policy, please contact
                      us at{" "}
                      <a href="mailto:privacy@chanzuckerberg.com">
                        privacy@chanzuckerberg.com
                      </a>
                      .
                      <p>
                        To comply with article 27 of the GDPR and the UK-GDPR,
                        we have appointed a representative who can accept
                        communications on behalf of CZIF and CZI LLC in relation
                        to personal data processing activities falling within
                        the scope of the GDPR or the UK-GDPR. If you wish to
                        contact them, their details are as follows:
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
                      </p>
                    </li>
                  </ol>
                </li>
              </ol>
            </main>
          </PrivacyStyle>
        </CommonStyle>
      </Layout>

      <Layout>
        <CommonStyle>
          <PrivacyStyle>
            {/* <SEO title="Terms of Service" /> */}
            <header>
              <CZILogo
                data-test-id="czi-logo"
                src={String(rawCZILogo)}
                alt="CZI logo"
              />
              <CellxgeneLogo
                data-test-id="cellxgene-logo"
                src={String(rawCellxgeneLogo)}
                alt="cellxgene logo"
              />
            </header>

            <h1 id="tos">Terms of Use</h1>

            <p>
              <span className="caps">
                ANY DISPUTE BETWEEN YOU AND US IS SUBJECT TO BINDING
                ARBITRATION, AS SET FORTH IN <a href="#section6">SECTION 6</a>.
                PLEASE READ THE ARBITRATION PROVISION IN THIS AGREEMENT AS IT
                AFFECTS YOUR RIGHTS UNDER THIS CONTRACT.
              </span>
            </p>

            <h2>Introduction</h2>

            <p>
              These Terms of Use (the “<strong>Terms</strong>”) apply to
              everyone who accesses and uses the version of cellxgene located at
              cellxgene.cziscience.com (the “<strong>Site</strong>”; those
              individuals who access and use it, “<strong>you</strong>”).
            </p>

            <p>
              The purpose of the Site is to provide a version of the interactive
              data explorer called cellxgene that enables fast visualizations of
              curated single-cell transcriptomics datasets (those datasets, the
              “<strong>Data</strong>”). By leveraging modern web development
              techniques to enable those visualizations, we hope to enable
              biologists and computational researchers to explore that Data.
            </p>

            <p>
              You agree that by accessing the Site, you are entering into a
              legally-binding contract with the Chan Zuckerberg Initiative
              Foundation, a 501(c)(3) nonprofit private foundation (
              <strong>“Provider,” “we,” “us”</strong>) and agree to be bound by
              these Terms as well as our{" "}
              <a href={ROUTES.PRIVACY}>Privacy Policy</a>. If you don’t agree to
              these Terms, or to our <a href={ROUTES.PRIVACY}>Privacy Policy</a>
              , don’t access the Site.
            </p>

            <p>
              Please read our full Terms and{" "}
              <a href={ROUTES.PRIVACY}>Privacy Policy</a> for complete details,
              but here is the key information you should know:
            </p>

            <ul>
              <li>
                This is a <strong>free service</strong> we provide in order to
                advance biomedical research.
              </li>
              <li>
                The data made available in the service is{" "}
                <strong>publicly available</strong> and is{" "}
                <strong>not personally identifiable</strong>.
              </li>
              <li>
                We collect <strong>basic analytics</strong> about your usage of
                the service (e.g. frequency, duration, IP address) so we know
                how it’s being used. This helps us improve the service as well
                as understand how it’s used.
              </li>
            </ul>

            <p className="caps">
              PLEASE READ THESE TERMS CAREFULLY AS USE OF THE SERVICE
              CONSTITUTES ACCEPTANCE OF THESE TERMS AND CONDITIONS.
            </p>

            <ol>
              <li>
                <h4>Use of the Services.</h4> You shall not otherwise access or
                use – or attempt to access or use – the Site to take any action
                that could harm us, the Site, or any third party, or use the
                Site in a manner that violates applicable law or violates the
                rights of others.
              </li>
              <li>
                <h4>Site Limitations.</h4> We may restrict or terminate your
                access to the Site at any time, without notice, and for any
                reason including for breach of these Terms. The Data on the Site
                has been compiled from a variety of sources, and is subject to
                change without notice. CZIF does not investigate, monitor or
                check Data for accuracy, appropriateness, completeness, or other
                reliability. As such, you agree that CZIF shall not be
                responsible for any Data or failure to include any Data or
                updates thereto. Your use of the Data is at your own risk.
              </li>

              <li>
                <h4>Disclaimer.</h4>{" "}
                <span className="caps">
                  THE SITE AND THE DATA ARE PROVIDED “AS IS” WITH ALL FAULTS,
                  AND WE AND OUR SERVICE PROVIDERS HEREBY DISCLAIM ALL
                  REPRESENTATIONS AND WARRANTIES, EXPRESS OR IMPLIED, WITH
                  RESPECT TO THE SITE AND THE DATA. THIS DISCLAIMER APPLIES TO
                  THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW.
                </span>
              </li>

              <li>
                <h4>Limitation of Liability.</h4>{" "}
                <span className="caps">
                  TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, CZIF AND
                  AFFILIATES (INCLUDING WITHOUT LIMITATION CHAN ZUCKERBERG
                  INITIATIVE, LLC; COLLECTIVELY, THE “PROTECTED PARTIES”) WILL
                  NOT BE LIABLE FOR ANY INDIRECT, INCIDENTAL, SPECIAL, PUNITIVE,
                  OR CONSEQUENTIAL DAMAGES OF ANY KIND (INCLUDING LOST PROFITS,
                  LOST DATA, BUSINESS INTERRUPTION, OR LOSS OF GOODWILL)
                  IRRESPECTIVE OF WHETHER SUCH DAMAGES ARISE FROM CLAIMS BROUGHT
                  IN CONTRACT, TORT, NEGLIGENCE, WARRANTY, STRICT LIABILITY, OR
                  ANY OTHER THEORY AT LAW OR IN EQUITY, AND EVEN IF ANY CZI LLC
                  PROTECTED PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
                  DAMAGES. WITHOUT LIMITING THE FOREGOING, TO THE MAXIMUM EXTENT
                  PERMITTED BY APPLICABLE LAW, IN NO EVENT WILL THE CZI LLC
                  PROTECTED PARTIES’ AGGREGATE LIABILITY ARISING OUT OF OR
                  RELATING TO THESE TERMS OR THE SERVICE EXCEED USD $100.
                </span>
              </li>

              <li>
                <h4>Indemnification.</h4>
                <ol className="section5">
                  <li>
                    You shall indemnify, defend and hold CZIF harmless from and
                    against, and shall pay all damages, costs, fees and expenses
                    (including reasonable attorneys’ fees and expenses) relating
                    to, any third party (including government entity) claim,
                    action, suit or other proceeding (a “Claim”) to the extent
                    arising from: (1) your gross negligence, willful misconduct
                    or fraud; and/or (2) any misrepresentation you make
                    regarding your permission to submit data to cellxgene for
                    public use.
                  </li>
                  <li>
                    Section 5.1 indemnification is conditioned upon, CZIF giving
                    you written notice of any such Claim, and giving you control
                    of the defense and settlement of any such Claim, and
                    cooperating with you in such defense. Notwithstanding
                    anything to the contrary, (1) CZIF may participate in
                    defense of such Claim with its own counsel at its own
                    expense and (2) you may not settle any Claim without CZIF’s
                    prior written consent, which will not be unreasonably
                    withheld, unless it unconditionally releases CZIF of all
                    liability, obligation, and fault.
                  </li>
                </ol>
              </li>

              <li id="section6">
                <h4>
                  Governing Law, Dispute Resolution and Mutual Agreement to
                  Arbitrate.
                </h4>
                <ol className="section6">
                  <li>
                    <h5>Final and Binding Arbitration.</h5> We endeavor and
                    trust that we will have a productive relationship but in the
                    unlikely event we have a dispute that we can’t resolve
                    between us, and it results in a legal dispute,{" "}
                    <strong>
                      <span className="caps">
                        BOTH YOU AND WE AGREE TO WAIVE OUR RESPECTIVE RIGHTS TO
                        RESOLUTION OF DISPUTES IN A COURT OF LAW BY A JUDGE OR
                        JURY AND AGREE TO RESOLVE ANY DISPUTE BY ARBITRATION,
                        WHICH WILL BE FINAL AND BINDING, AS SET FORTH BELOW.
                      </span>
                    </strong>
                  </li>
                  <li>
                    <h5>Dispute Resolution.</h5> In the unlikely event we have a
                    dispute arising out of or related to the use of the Site
                    (“Dispute”) that we can’t resolve between us, you and we
                    agree that we shall (in good faith) meet and attempt to
                    resolve the Dispute within thirty (30) days. If the Dispute
                    is not resolved during such time period, then you and a
                    representative of CZIF shall (in good faith) meet and
                    attempt to resolve the Dispute through non-binding mediation
                    with a mutually agreed upon mediator within thirty (30)
                    additional days.
                  </li>
                  <li>
                    <h5>Mutual Agreement to Arbitrate.</h5> If the Dispute is
                    not resolved within such time period, the Dispute shall be
                    resolved per the following arbitration terms. As the
                    exclusive, final and binding means of initiating adversarial
                    proceedings, you agree that it be resolved fully and finally
                    by neutral and binding arbitration administered by JAMS in
                    San Mateo County, California, in accordance with its
                    Streamlined Arbitration Rules & Procedures, the Federal
                    Arbitration Act, and the substantive laws of the State of
                    California, exclusive of conflict or choice of law rules.
                    In-person proceedings will take place in San Mateo County,
                    California and your reasonable and documented travel
                    expenses will be paid by CZIF. The arbitrator shall have the
                    power to award any type of relief that would be available in
                    a court of competent jurisdiction and will issue a written
                    decision at the end of the arbitration, which will be final
                    and binding. Judgment on any award rendered in any such
                    arbitration may be entered in any court having jurisdiction
                    in San Mateo County, California.
                  </li>
                  <li>
                    <h5>Where inapplicable, Choice of Law and Venue.</h5> This
                    Agreement and any Disputes will be governed, controlled, and
                    interpreted by and under the laws of the State of
                    California, without giving effect to any conflicts of laws
                    principles that require the application of the law of a
                    different state. Notwithstanding the foregoing, to the
                    extent such laws are inconsistent with the Federal
                    Arbitration Act, the Federal Arbitration Act will govern.
                    Any dispute that is not subject to arbitration (e.g., if
                    arbitration is deemed unenforceable or inapplicable) shall
                    be, and any judgement on any arbitration award may be,
                    brought in the U.S. District Court for the Northern District
                    of California or a state court located in San Mateo County,
                    California.
                  </li>
                </ol>
              </li>

              <li id="section7">
                <h4>General Terms.</h4>
                <ol className="section7">
                  <li>
                    If any provision in these Terms is held invalid or
                    unenforceable, the other provisions will remain enforceable,
                    and the invalid or unenforceable provision will be modified
                    to a valid and enforceable provision that most accurately
                    reflects the parties intentions.
                  </li>
                  <li>
                    Any waiver or failure to enforce any of these Terms on one
                    occasion will not be deemed a waiver of any other provision
                    or of that provision on any other occasion.
                  </li>
                  <li>
                    You may not assign or transfer any rights or obligations
                    under these Terms without our consent. However, you agree
                    that we may assign these Terms in connection with a
                    reorganization, or to a successor or assign that agrees to
                    assume our obligations under these Terms (and{" "}
                    <a href={ROUTES.PRIVACY}>Privacy Policy</a>) without your
                    consent.
                  </li>
                </ol>
              </li>

              <li>
                <h4>How To Contact Us.</h4> Notice under these Terms must be in
                writing and deemed to have been given on the date delivered by a
                nationally recognized express mail service, such as Federal
                Express, or by certified and registered mail (signature for
                receipt required) to CZIF as follows:
                <br />
                <br />
                <address>
                  Chan Zuckerberg Initiative Foundation <br />
                  314 Lytton Ave.
                  <br />
                  Palo Alto, CA 94301
                  <br />
                  Attn: General Counsel <br />
                  Email: courtesy copy:{" "}
                  <a href="mailto:legalczi1@chanzuckerberg.com">
                    legalczi1@chanzuckerberg.com
                  </a>{" "}
                  (email does not constitute notice)
                </address>
              </li>

              <li>
                <h4>Data Submission.</h4> If you would like to submit Data to
                the Site, please contact{" "}
                <a href="mailto:cellxgene@chanzuckerberg.com">
                  cellxgene@chanzuckerberg.com
                </a>
                . Keep in mind that – if the Data you wish to submit is not
                already published in a public archive such as GEO, as an open
                link in a publication, or in a github repository – you will be
                required to grant CZIF permission to use, display and create
                derivative works (e.g. visualizations) of the Data for purposes
                of offering the Site, and must therefore have the authority to
                give that permission.
              </li>
            </ol>
          </PrivacyStyle>
        </CommonStyle>
      </Layout>
    </Layout>
  );
};

export default PreviewPolicies;
