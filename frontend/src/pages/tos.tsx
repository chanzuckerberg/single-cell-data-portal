import Head from "next/head";
import Image from "next/image";
import { ROUTES } from "src/common/constants/routes";
import rawCellxgeneLogo from "src/components/common/staticPages/cellxgene.png";
import {
  CommonStyle,
  Layout,
  TOSStyle,
} from "src/components/common/staticPages/style";

const ToS = (): JSX.Element => {
  return (
    <Layout>
      <CommonStyle>
        <TOSStyle>
          <Head>
            <title>Terms of Service - CZ CELLxGENE</title>
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

          <h1 id="tos">Terms of Use</h1>

          <p>Last Updated: December 15, 2022.</p>

          <p>
            <span className="caps">
              ANY DISPUTE BETWEEN YOU AND US IS SUBJECT TO BINDING ARBITRATION,
              AS SET FORTH IN <a href="#section6">SECTION 6</a>. PLEASE READ THE
              ARBITRATION PROVISION IN THIS AGREEMENT AS IT AFFECTS YOUR RIGHTS
              UNDER THIS CONTRACT.
            </span>
          </p>

          <br />

          <h2>Introduction</h2>

          <p>
            These Terms of Use (the “<strong>Terms</strong>”) apply to everyone
            who accesses and uses the version of CELLxGENE located at
            cellxgene.cziscience.com (the “<strong>Site</strong>”; those
            individuals who access and use it, “<strong>you</strong>”).
          </p>

          <p>
            The purpose of the Site is to provide a version of the interactive
            data explorer called CELLxGENE that enables fast visualizations of
            curated single-cell transcriptomics datasets (those datasets, the “
            <strong>Data</strong>”). By leveraging modern web development
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
            these Terms, or to our <a href={ROUTES.PRIVACY}>Privacy Policy</a>,
            don’t access the Site.
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
              We use the privacy-friendly Plausible service to collect{" "}
              <strong>basic analytics</strong> about Site traffic (e.g. the
              number of visitors and page views) so we know how it’s being used.
              This helps us improve the service as well as understand how it’s
              used.
            </li>
          </ul>

          <p className="caps">
            PLEASE READ THESE TERMS CAREFULLY AS USE OF THE SERVICE CONSTITUTES
            ACCEPTANCE OF THESE TERMS AND CONDITIONS.
          </p>

          <ol>
            <li>
              <h4>Use of the Services.</h4> You shall not otherwise access or
              use – or attempt to access or use – the Site to take any action
              that could harm us, the Site, or any third party, or use the Site
              in a manner that violates applicable law or violates the rights of
              others.
            </li>
            <li>
              <h4>Site Limitations.</h4> We may restrict or terminate your
              access to the Site at any time, without notice, and for any reason
              including for breach of these Terms. The Data on the Site has been
              compiled from a variety of sources, and is subject to change
              without notice. CZIF does not investigate, monitor or check Data
              for accuracy, appropriateness, completeness, or other reliability.
              As such, you agree that CZIF shall not be responsible for any Data
              or failure to include any Data or updates thereto. Your use of the
              Data is at your own risk.
            </li>

            <li>
              <h4>Disclaimer.</h4>{" "}
              <span className="caps">
                THE SITE AND THE DATA ARE PROVIDED “AS IS” WITH ALL FAULTS, AND
                WE AND OUR SERVICE PROVIDERS HEREBY DISCLAIM ALL REPRESENTATIONS
                AND WARRANTIES, EXPRESS OR IMPLIED, WITH RESPECT TO THE SITE AND
                THE DATA. THIS DISCLAIMER APPLIES TO THE MAXIMUM EXTENT
                PERMITTED BY APPLICABLE LAW.
              </span>
            </li>

            <li>
              <h4>Limitation of Liability.</h4>{" "}
              <span className="caps">
                TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, CZIF AND
                AFFILIATES (INCLUDING WITHOUT LIMITATION CHAN ZUCKERBERG
                INITIATIVE, LLC; COLLECTIVELY, THE “PROTECTED PARTIES”) WILL NOT
                BE LIABLE FOR ANY INDIRECT, INCIDENTAL, SPECIAL, PUNITIVE, OR
                CONSEQUENTIAL DAMAGES OF ANY KIND (INCLUDING LOST PROFITS, LOST
                DATA, BUSINESS INTERRUPTION, OR LOSS OF GOODWILL) IRRESPECTIVE
                OF WHETHER SUCH DAMAGES ARISE FROM CLAIMS BROUGHT IN CONTRACT,
                TORT, NEGLIGENCE, WARRANTY, STRICT LIABILITY, OR ANY OTHER
                THEORY AT LAW OR IN EQUITY, AND EVEN IF ANY CZI LLC PROTECTED
                PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
                WITHOUT LIMITING THE FOREGOING, TO THE MAXIMUM EXTENT PERMITTED
                BY APPLICABLE LAW, IN NO EVENT WILL THE CZI LLC PROTECTED
                PARTIES’ AGGREGATE LIABILITY ARISING OUT OF OR RELATING TO THESE
                TERMS OR THE SERVICE EXCEED USD $100.
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
                  arising from: (1) your gross negligence, willful misconduct or
                  fraud; and/or (2) any misrepresentation you make regarding
                  your permission to submit data to CELLxGENE for public use.
                </li>
                <li>
                  Section 5.1 indemnification is conditioned upon, CZIF giving
                  you written notice of any such Claim, and giving you control
                  of the defense and settlement of any such Claim, and
                  cooperating with you in such defense. Notwithstanding anything
                  to the contrary, (1) CZIF may participate in defense of such
                  Claim with its own counsel at its own expense and (2) you may
                  not settle any Claim without CZIF’s prior written consent,
                  which will not be unreasonably withheld, unless it
                  unconditionally releases CZIF of all liability, obligation,
                  and fault.
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
                  <h5>Final and Binding Arbitration.</h5> We endeavor and trust
                  that we will have a productive relationship but in the
                  unlikely event we have a dispute that we can’t resolve between
                  us, and it results in a legal dispute,{" "}
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
                  (“Dispute”) that we can’t resolve between us, you and we agree
                  that we shall (in good faith) meet and attempt to resolve the
                  Dispute within thirty (30) days. If the Dispute is not
                  resolved during such time period, then you and a
                  representative of CZIF shall (in good faith) meet and attempt
                  to resolve the Dispute through non-binding mediation with a
                  mutually agreed upon mediator within thirty (30) additional
                  days.
                </li>
                <li>
                  <h5>Mutual Agreement to Arbitrate.</h5> If the Dispute is not
                  resolved within such time period, the Dispute shall be
                  resolved per the following arbitration terms. As the
                  exclusive, final and binding means of initiating adversarial
                  proceedings, you agree that it be resolved fully and finally
                  by neutral and binding arbitration administered by JAMS in San
                  Mateo County, California, in accordance with its Streamlined
                  Arbitration Rules & Procedures, the Federal Arbitration Act,
                  and the substantive laws of the State of California, exclusive
                  of conflict or choice of law rules. In-person proceedings will
                  take place in San Mateo County, California and your reasonable
                  and documented travel expenses will be paid by CZIF. The
                  arbitrator shall have the power to award any type of relief
                  that would be available in a court of competent jurisdiction
                  and will issue a written decision at the end of the
                  arbitration, which will be final and binding. Judgment on any
                  award rendered in any such arbitration may be entered in any
                  court having jurisdiction in San Mateo County, California.
                </li>
                <li>
                  <h5>Where inapplicable, Choice of Law and Venue.</h5> This
                  Agreement and any Disputes will be governed, controlled, and
                  interpreted by and under the laws of the State of California,
                  without giving effect to any conflicts of laws principles that
                  require the application of the law of a different state.
                  Notwithstanding the foregoing, to the extent such laws are
                  inconsistent with the Federal Arbitration Act, the Federal
                  Arbitration Act will govern. Any dispute that is not subject
                  to arbitration (e.g., if arbitration is deemed unenforceable
                  or inapplicable) shall be, and any judgement on any
                  arbitration award may be, brought in the U.S. District Court
                  for the Northern District of California or a state court
                  located in San Mateo County, California.
                </li>
              </ol>
            </li>

            <li id="section7">
              <h4>General Terms.</h4>
              <ol className="section7">
                <li>
                  If any provision in these Terms is held invalid or
                  unenforceable, the other provisions will remain enforceable,
                  and the invalid or unenforceable provision will be modified to
                  a valid and enforceable provision that most accurately
                  reflects the parties intentions.
                </li>
                <li>
                  Any waiver or failure to enforce any of these Terms on one
                  occasion will not be deemed a waiver of any other provision or
                  of that provision on any other occasion.
                </li>
                <li>
                  You may not assign or transfer any rights or obligations under
                  these Terms without our consent. However, you agree that we
                  may assign these Terms in connection with a reorganization, or
                  to a successor or assign that agrees to assume our obligations
                  under these Terms (and{" "}
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
                801 Jefferson Ave.
                <br />
                Redwood City, CA 94063
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
              <h4>Data Submission.</h4> If you would like to submit Data to the
              Site, please contact{" "}
              <a href="mailto:cellxgene@chanzuckerberg.com">
                cellxgene@chanzuckerberg.com
              </a>
              . Keep in mind that – if the Data you wish to submit is not
              already published in a public archive such as GEO, as an open link
              in a publication, or in a github repository – you will be required
              to grant CZIF permission to use, display and create derivative
              works (e.g. visualizations) of the Data for purposes of offering
              the Site, and must therefore have the authority to give that
              permission.
            </li>
          </ol>
        </TOSStyle>
      </CommonStyle>
    </Layout>
  );
};

export default ToS;
